import CifFile as cf
import sqlite3
import logging
import re
import pymatgen.core as pmg
import numpy as np

## SHARED METHODS IN CORE CLASS
class core:

    ION_REGEX = re.compile("([A-Za-z]{,2})(\d*)(\+|-)")
    RESULT_CONVERTER = lambda self, x: "1" if x == "" else x

    def interpretIon(self, ion):
            """
                Interprets an ion in string format. The ion must be of format Symbol-OS Digit-Sign, otherwise the interpretation will not work.

                Returns a tuple of (element, oxidation_state)
            """
            result = self.ION_REGEX.search(ion)
            if result is None:
                raise Exception(f"Cannot interpret ion ({ion}) sucessfully.")
            else:
                return (result[1], int(result[3] + self.RESULT_CONVERTER(result[2])))

## DATABASE CLASS

class BVDatabase:
    """
        A class representing a connection to a bond valence parameter database. Contains all methods required to communitate with the database.
    """

    CORE = core()
    ION_REGEX = re.compile("([A-Za-z]{,2})(\d*)(\+|-)")
    RESULT_CONVERTER = lambda self, x: "1" if x == "" else x

    def __init__(self, dbLocation:str):
        """
            Initialise the connection to the database using its location.
        """
        self.conn = sqlite3.connect(dbLocation)
        self.cursor = self.conn.cursor()

    def commit(self):
        """
            Commit changes to the database.
        """
        self.conn.commit()
        
    def close(self):
        """
            Close the database, first commiting any remaing changes.
        """
        self.commit()
        self.conn.close()

    def execute(self, command:str, arguments = None):
        """
            Execute an SQL command to the database. Can optionally take arguments.
        """
        if arguments is None:
            self.cursor.execute(command)
        else:
            self.cursor.execute(command, arguments)

    def fetchall(self):
        """
            Method for fetching all results from the database.
        """
        return self.cursor.fetchall()

    def fetchone(self):
        """
            Method for fetching one result from the database.
        """
        return self.cursor.fetchone()
    
    def extractFetchone(self):
        """
            Method for removing the tuple around the results from fetchone().
        """
        result = self.fetchone()
        if result is not None:
            return result[0]
        else:
            return None
    
    def lastRowId(self):
        """
            Method for getting the id of the last row changed on the database.
        """
        return self.cursor.lastrowid

    def createDatabase(self):
        """
            Executes a SQL script to reset the database.
        """

        self.cursor.executescript('''
            CREATE TABLE Ion (
                id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
                symbol VARCHAR(2),
                os INTEGER(1),
                radii REAL,
                softness REAL                     
            );

            CREATE TABLE BVParam (
                id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
                ion1 INTEGER,
                ion2 INTEGER,
                r0 REAL,
                b REAL,
                ib REAL,
                cn REAL,
                r_cutoff REAL
            );
            ''')
        
    def resetDatabase(self):
        """
            Method that resets the database by dropping all tables and executing a SQL script to create the database again.
        """

        # Checks the user really wants to reset the database
        print("Are you sure you want to reset the database? (No data can be retrieved)")
        ans = input().lower().strip()
        if ans != "y" and ans != 'yes':
            return
      
        # Code for dropping tables lifted from db-reset.py
        self.execute("SELECT 'drop table ' || name || ';' FROM sqlite_master WHERE type = 'table'")
        commands = self.fetchall()

        for (command,) in commands:
            if "sqlite_sequence" in command:
                logging.info("Skipped command {c} (Intentional)".format(c=command))
                continue
            self.execute(command)
            logging.info("Succesfully Executed {c}".format(c=command))

        # Creates the new database
        self.createDatabase()

    def getOrInsertIon(self, ion:str, os:int):
        """
            Method that tries to find an ion in the database. If it is not found, a new database entry is added for that ion. In any case, the id of the ion's entry is returned.
        """

        self.execute("SELECT id FROM Ion WHERE symbol = ? AND os = ?", (ion, os, ))

        ans = self.extractFetchone()
        if ans == None:
            self.execute("INSERT OR IGNORE INTO Ion (symbol, os) VALUES (?,?)", (ion, os, ))
            return self.lastRowId()
        else:
            return ans

    def createEntry(self, ion1:str, os1:int, ion2:str, os2:int, r0:float, b:float):
        """
            Creates a database entry for a BV parameter set. Only basic information, included in Brown's cif files is set in this method.
        """

        ion1Id = self.getOrInsertIon(ion1, os1)
        ion2Id = self.getOrInsertIon(ion2, os2)

        self.execute("SELECT id FROM BVParam WHERE ion1=? AND ion2=?", (ion1Id, ion2Id))
        if self.extractFetchone() == None:
            ib = 1/b
            self.execute("INSERT INTO BVParam (ion1, ion2, r0, b, ib) VALUES (?,?,?,?,?)", (ion1Id, ion2Id, r0, b, ib))

        return self.lastRowId()
    
    def addInfo(self, paramId:int, cn:float, rCutoff:float):
        """
            Updates the database entry to add information on the coordination number and radius cutoff, information which is only included in the softBV parameters.
        """

        self.execute("UPDATE BVParam SET cn = ?, r_cutoff = ? WHERE id = ?", (cn, rCutoff, paramId, ))

    
    def getParams(self, ion1:str, ion2:str):
        """
            Finds the parameters for a combination of two ions and returns them. The input ions should be in the format of Symbol-OS Digit-Sign.

            Returns a tuple of (r0, ib)
        """
        ion1, ion2 = self.CORE.interpretIon(ion1), self.CORE.interpretIon(ion2)
        self.execute("SELECT r0, ib FROM BVParam JOIN Ion i1 JOIN Ion i2 On BVParam.ion1 = i1.id AND BVParam.ion2 = i2.id WHERE i1.symbol = ? AND i1.os = ? AND i2.symbol = ? AND i2.os = ?", (ion1[0], ion1[1], ion2[0], ion2[1],))

        result = self.fetchall()
        if result is None or len(result) == 0:
            raise Exception(f"The combination of ions ({ion1}, {ion2}) are not on the BV Parameters Database")
        elif len(result) != 1:
            logging.warn(f"Multiple different database entries for the same ions ({ion1}, {ion2})")
        
        return result[0]


def readCif(fileLocation:str):
    return cf.ReadCif(fileLocation).first_block()

def fileToDb(fileIn:str, fileOut:str):
    """
        Converts any recognised file into a bond valence parameter database
    """
    if fileIn[-4:] == ".dat":
        datToDb(fileIn, fileOut)
    elif fileIn[-4:] == ".cif":
        cifToDb(fileIn, fileOut) 
    else:
        print("Data format unrecognised - Can only read .cif and .dat files")

def cifToDb(fileIn:str, fileOut:str):
    """
        Convert the bond valence cif file provied on the ICuR website into a local database format for easy access.
    """
    
    cif = readCif(fileIn)
    bvTable = cif.GetLoop("_valence_param_atom_1")
    db = BVDatabase(fileOut)

    for row in bvTable:
        db.createEntry(row[0], int(row[1]), row[2], int(row[3]), float(row[4]), float(row[5]))

    db.close()

def datToDb(fileIn:str, fileOut:str):
    """
        Convert the bond valence dat file from softBV into a local database format
    """

    db = BVDatabase(fileOut)
    db.resetDatabase()

    with open(fileIn, "r") as data:
        entries = data.readlines()
        
        started = False
        for entry in entries:
            if started and entry == "DATA_END\n":
                break
            elif started:
                row =  entry.split()
                paramId = db.createEntry(row[0], int(row[1]), row[2], int(row[3]), float(row[4]), float(row[5]))
                db.addInfo(paramId, float(row[6]), float(row[7]))
            elif entry == "DATA_START\n":
                started = True
            else:
                continue

    db.close()

def createInputFromCif(fileIn:str, fileOut:str, conductor:str):
    
    CORE = core()

    # Create a pymatgen object
    struct = pmg.Structure.from_file(fileIn)
    conductorSymb, conductorOS = CORE.interpretIon(conductor)

    # Open output file
    with open(fileOut, "w") as f:

        # Give information on conductor choice and the lattice
        f.write(f"{conductorSymb}\t{conductorOS}\n")
        f.write(f"{struct.lattice.a}\t{struct.lattice.b}\t{struct.lattice.c}\t{struct.lattice.alpha}\t{struct.lattice.beta}\t{struct.lattice.gamma}\n")
        f.write(f"{struct.lattice.volume}\n")
        for i in range(3):
            for j in range (3):
                f.write(f"{struct.lattice.matrix[i][j]}\t")
            f.write("\n")

        # For every site, add label, element, os and cartesian coords
        for site in struct.sites:

            f.write(site.label + "\t")

            if isinstance(site.species, pmg.Composition):
                elements = site.species.elements
                if len(elements) != 1:
                    f.write("##DISORDERED SITE - AMEND MANUALLY##\t")
                else:
                    f.write(f"{elements[0].element}\t{elements[0].oxi_state}\t")
            
            elif isinstance(site.species, pmg.Species):
                f.write(f"{site.species.symbol}\t{site.species.oxidation_state}\t")
            else:
                f.write("##DISORDERED SITE - AMEND MANUALLY##\t")
            
            f.write(f"{site.coords[0]}\t{site.coords[1]}\t{site.coords[2]}\n")


    
createInputFromCif("cif-files/Ternary Fluorides/EntryWithCollCode152949 (PbSnF4).cif", "pbsnf4.inp", "F-")

# datToDb("cif-files/database_binary.dat", "soft-bv-params.sqlite3")
# db = BVDatabase("soft-bv-params.sqlite3")
# print(db.getParams("Pb2+","O2-"))
        