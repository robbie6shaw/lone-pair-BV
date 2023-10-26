import CifFile as cf
import pandas as pd
import numpy as np
import sqlite3
import logging
import re
import pymatgen.core as pmg
import numpy as np
import collections

## SHARED METHODS IN CORE CLASS
class core:

    ION_REGEX = re.compile("([A-Za-z]{,2})(\d*)(\+|-)")
    RESULT_CONVERTER = lambda self, x: "1" if x == "" else x
    LONE_PAIR_ELEMENTS = ["Pb", "Sn", "Bi", "Sb", "Tl"]
    ion = collections.namedtuple("Ion", ["element", "os"])
    bvparam = collections.namedtuple("BVParam", ['r0', 'ib', 'cn', 'r_cutoff', 'i1r', 'i2r', 'rmin', 'd0'])

    # def interpretIon(self, ion):
    #         """
    #             Interprets an ion in string format. The ion must be of format Symbol-OS Digit-Sign, otherwise the interpretation will not work.

    #             Returns a tuple of (element, oxidation_state)
    #         """
    #         result = self.ION_REGEX.search(ion)
    #         if result is None:
    #             raise Exception(f"Cannot interpret ion ({ion}) sucessfully.")
    #         else:
    #             return self.ion(element = result[1], os = int(result[3] + self.RESULT_CONVERTER(result[2])))
    
    # def ion_to_str(self, ion:tuple|ion):
    #     """
    #         Takes an ion in a tuple format and returns it in string format.
    #     """
    #     os = int(ion[1])
    #     return ion[0] + (str(abs(os)) if abs(os) > 1 else "") + ("+" if os > 0 else "-")
            
    def hasLonePair(self, element:str):

        return element in self.LONE_PAIR_ELEMENTS
            

class Ion:

    ION_REGEX = re.compile("([A-Za-z]{,2})(\d*)(\+|-)")
    LONE_PAIR_ELEMENTS = ["Pb", "Sn", "Bi", "Sb", "Tl"]
    one_reducer = lambda self, x: "1" if x == "" else x

    def __init__(self, element:str, ox_state:int):

        if isinstance(element, str):
            self.element = element
        else:
            raise TypeError("The element field of an Ion should be a string")
        
        self.ox_state = int(ox_state)
        self.string = self.element + (str(abs(self.ox_state)) if abs(self.ox_state) > 1 else "") + ("+" if self.ox_state > 0 else "-")
        self.radius = None

    def __str__(self):
        return self.string
    
    def __eq__(self, other):
        if isinstance(other, Ion):
            return self.element == other.element and self.ox_state == other.ox_state
        elif isinstance(other, str):
            other = Ion.from_string(other)
            return self.element == other.element and self.ox_state == other.ox_state
        else:
            return False
        
    def __hash__(self):
        return hash((self.element, self.ox_state))

    @classmethod
    def from_string(cls, ion_string:str):
        result = cls.ION_REGEX.search(ion_string)
        
        if result is None:
            raise Exception(f"Cannot interpret the ion string sucessfully - {ion_string}")
        else:
            return Ion(element = result[1], os = int(result[3] + cls.one_reducer(result[2])))
        
    def possible_lone_pair(self):
        return self.element in self.LONE_PAIR_ELEMENTS

## DATABASE CLASS

class BVDatabase:
    """
        A class representing a connection to a bond valence parameter database. Contains all methods required to communitate with the database.
    """

    CORE = core()
    BVParams = collections.namedtuple("BVParams", ['r0', 'b', 'ib', 'cn', 'r_cutoff', 'i1_radius', 'i2_radius', 'i1_softness', 'i2_softness', 'i1_period', 'i2_period','i1_block', 'i2_block'])

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

    def fetch_all(self):
        """
            Method for fetching all results from the database.
        """
        return self.cursor.fetchall()
    
    def fetch_one(self):
        """
            Method for removing the tuple around the results from fetchone().
        """
        result = self.cursor.fetchone()
        if result is not None:
            return result[0]
        else:
            return None
    
    def last_row_id(self):
        """
            Method for getting the id of the last row changed on the database.
        """
        return self.cursor.lastrowid

    def create_database(self):
        """
            Executes a SQL script to reset the database.
        """

        self.cursor.executescript('''
            CREATE TABLE Ion (
                id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
                symbol VARCHAR(2),
                atomic_no INT,
                os INTEGER(1),
                radii REAL,
                softness REAL,
                period INT,
                p_group INT,
                block INT                
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
        
    def reset_database(self):
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
        commands = self.fetch_all()

        for (command,) in commands:
            if "sqlite_sequence" in command:
                logging.info("Skipped command {c} (Intentional)".format(c=command))
                continue
            self.execute(command)
            logging.info("Succesfully Executed {c}".format(c=command))

        # Creates the new database
        self.create_database()

    def get_or_insert_ion(self, ion:Ion):
        """
            Method that tries to find an ion in the database. If it is not found, a new database entry is added for that ion. In any case, the id of the ion's entry is returned.
        """

        self.execute("SELECT id FROM Ion WHERE symbol = ? AND os = ?", (ion.element, ion.ox_state, ))

        ans = self.fetch_one()
        if ans == None:
            self.execute("INSERT OR IGNORE INTO Ion (symbol, os) VALUES (?,?)", (ion.element, ion.ox_state, ))
            return self.last_row_id()
        else:
            return ans

    def create_entry(self, ion1:Ion, ion2:Ion, r0:float, b:float):
        """
            Creates a database entry for a BV parameter set. Only basic information, included in Brown's cif files is set in this method.
        """

        ion1Id = self.get_or_insert_ion(ion1)
        ion2Id = self.get_or_insert_ion(ion2)

        self.execute("SELECT id FROM BVParam WHERE ion1=? AND ion2=?", (ion1Id, ion2Id))
        if self.fetch_one() == None:
            ib = 1/b
            self.execute("INSERT INTO BVParam (ion1, ion2, r0, b, ib) VALUES (?,?,?,?,?)", (ion1Id, ion2Id, r0, b, ib))

        return self.last_row_id()
    
    def update_bv_info(self, paramId:int, cn:float, rCutoff:float):
        """
            Updates the bond valence parameters database entry to add information on the coordination number and radius cutoff, information which is only included in the softBV parameters.
        """

        self.execute("UPDATE BVParam SET cn = ?, r_cutoff = ? WHERE id = ?", (cn, rCutoff, paramId, ))

    def update_ion_info(self, ionId:int, radii:float, softness:float, period:int, group:int, block:int, atomicNo:int):
        """
            Updates the ion database entry to add information on the radius, softness and principle quantum number, only included in the softBV parameters.
        """

        self.execute("UPDATE Ion SET radii = ?, softness = ?, period = ?, p_group = ?, block = ?, atomic_no = ? WHERE id = ?", (radii, softness, period, group, block, atomicNo, ionId))

    def rmin(self, softness1:float, softness2:float, r0:float, b:float, osCation:int, cn:float):
        x = (0.9185 + 0.2285 * abs(softness1 - softness2))*r0
        y = b*np.log(abs(osCation)/cn)
        return x - y
    
    def d0(self, b:float, os1:int, os2:int, block1:int, rmin:float, period1:int, period2:int):
        if block1 <= 1: c = 1
        else: c = 2 

        return ((b**2)/2 * 14.4 * (c*(abs(os1*os2))**(1/c)) / (rmin*np.sqrt(period1 * period2)))
    
    def _params_error_check(self, ion1:Ion, ion2:Ion):
        result = self.fetch_all()
        if result is None or len(result) == 0:
            raise Exception(f"The combination of ions ({ion1}, {ion2}) are not on the BV Parameters Database")
        elif len(result) != 1:
            logging.warn(f"Multiple different database entries for the same ions ({ion1}, {ion2})")
        return result[0]

    def get_bv_params(self, ion1:Ion, ion2:Ion, bvse:bool = False):
        """
            Finds the parameters for a combination of two ions and returns them. The input ions should be in the format of Symbol-OS Digit-Sign.

            Returns a tuple of (r0, ib)
        """
        
        if (ion1.ox_state * ion2.ox_state) < 0:
            #                    0   1  2   3   4         5         6         7            8            9          10         11        12
            self.execute("SELECT r0, b, ib, cn, r_cutoff, i1.radii, i2.radii, i1.softness, i2.softness, i1.period, i2.period, i1.block, i2.block FROM BVParam JOIN Ion i1 JOIN Ion i2 On BVParam.ion1 = i1.id AND BVParam.ion2 = i2.id WHERE (i1.symbol = ? AND i1.os = ? AND i2.symbol = ? AND i2.os = ?) OR (i2.symbol = ? AND i2.os = ? AND i1.symbol = ? AND i1.os = ?)", (ion1.element, ion1.ox_state, ion2.element, ion2.ox_state, ion1.element, ion1.ox_state, ion2.element, ion2.ox_state,))
            
            result = self._params_error_check(ion1, ion2)

            if bvse:
                if ion1.ox_state > 0: cationOs = ion1.ox_state
                else: cationOs = ion2.ox_state
                rmin = self.rmin(result[7], result[8], result[0], result[1], cationOs, result[3])
                d0 = self.d0(result[1], ion1.ox_state, ion2.ox_state, result[11], rmin, result[9], result[10])
                return core.bvparam(r0 = result[0], ib = result[2], cn = result[3], r_cutoff = result[4], i1r = result[5], i2r = result[6], rmin = rmin, d0 = d0)
            else:
                return core.bvparam(r0 = result[0], ib = result[2], cn = result[3], r_cutoff = result[4], i1r = None, i2r = None, rmin = None, d0 = None)
            
        elif bvse:
            self.execute("SELECT i1.radii, i2.radii FROM Ion i1 JOIN Ion i2 WHERE (i1.symbol = ? AND i1.os = ?) OR (i2.symbol = ? AND i2.ion = ?)", ion1[0], int(ion1[1]), ion2[0], int(ion2[2]))

            result = self._params_error_check(ion1, ion2)

            return core.bvparam(r0=None, ib=None, cn=None, r_cutoff=None, i1r=result[0], i2r=result[1], rmin=None, d0=None)
        else:
            return None

        
    def get_atomic_no(self, element:str):
        """
            Function to find the atomic number of a particular element, given the symbol.
        """

        self.execute("SELECT atomic_no FROM Ion WHERE symbol = ?", (element,))
        return self.fetch_one()
    
    def get_radius(self, ion:Ion):
        """
            Function to find the ionic radius of a particular element, given the symbol
        """
        self.execute("SELECT radii FROM Ion WHERE symbol = ? AND os = ?", (ion.element, ion.ox_state, ))
        return self.fetch_one()
    
    def get_period(self, ion:Ion):
        """
            Function to find the period of a particular element, given the symbol
        """
        self.execute("SELECT period FROM Ion WHERE symbol = ? AND os = ?", (ion.element, ion.ox_state, ))
        return self.fetch_one()


def readCif(fileLocation:str):
    return cf.ReadCif(fileLocation).first_block()

def fileToDb(fileIn:str, fileOut:str):
    """
        Converts any recognised file into a bond valence parameter database
    """
    if fileIn[-4:] == ".dat":
        bvDatToDb(fileIn, fileOut)
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
        db.create_entry(row[0], int(row[1]), row[2], int(row[3]), float(row[4]), float(row[5]))

    db.close()

def bvDatToDb(fileIn:str, fileOut:str):
    """
        Convert the bond valence dat file from softBV into a local database format
    """

    db = BVDatabase(fileOut)
    db.reset_database()

    with open(fileIn, "r") as data:
        entries = data.readlines()
        
        started = False
        for entry in entries:
            if started and entry == "DATA_END\n":
                break
            elif started:
                row =  entry.split()
                paramId = db.create_entry(row[0], int(row[1]), row[2], int(row[3]), float(row[4]), float(row[5]))
                db.update_bv_info(paramId, float(row[6]), float(row[7]))
            elif entry == "DATA_START\n":
                started = True
            else:
                continue

    db.close()

def ionDatToDb(fileIn:str, fileOut:str):
    """
        Convert ion database file from softBV into a local database format
    """

    db = BVDatabase(fileOut)
    
    with open(fileIn, "r") as data:
        entries = data.readlines()

        started = False
        for entry in entries:
            if started and entry == "DATA_END\n":
                break
            elif started:
                row = entry.split()
                ionId = db.get_or_insert_ion(Ion(row[1], int(row[2])))
                db.update_ion_info(ionId, float(row[8]), float(row[9]), int(row[5]), int(row[6]), int(row[7]), int(row[0]))
            elif entry == "DATA_START\n":
                started = True
            else:
                continue

    db.close()


def create_input_from_cif(fileIn:str, fileOut:str, conductor:str):
    
    CORE = core()

    # Create a pymatgen object
    struct = pmg.Structure.from_file(fileIn)
    conductor = Ion.from_string(conductor)

    if fileIn[-4:] != ".cif":
        logging.error(f"Incorrect file format. The input file should be a cif file and not a {fileIn[-3:]} file")
    elif fileOut[-4:] != ".inp":
        logging.error(f"Incorrect file format. The output file should be a inp file and not a {fileIn[-3:]} file")

    # Open output file
    with open(fileOut, "w") as f:

        # Give information on conductor choice and the lattice
        f.write(f"{conductor.element}\t{conductor.ox_state}\n")
        f.write(f"{struct.lattice.a}\t{struct.lattice.b}\t{struct.lattice.c}\t{struct.lattice.alpha}\t{struct.lattice.beta}\t{struct.lattice.gamma}\n")
        f.write(f"{struct.lattice.volume}\n")
        for i in range(3):
            for j in range (3):
                f.write(f"{struct.lattice.matrix[i][j]}\t")
            f.write("\n")

        f.write("sym_label\tp1_label\telement\tos\tlp\ta\tb\tc\n")
        siteDict = {}
        osWarning = False

        # Highly inefficient, highly simple
        for site in struct.sites:
            siteDict[site.label] = 0

        # For every site, add label, element, os and cartesian coords
        for site in struct.sites:
            f.write(f"{site.label}\t{site.label}.{siteDict[site.label]}\t")
            siteDict[site.label] += 1
            if isinstance(site.species, pmg.Composition):
                
                elements = site.species.elements
                if len(elements) != 1:
                    f.write("##DISORDERED SITE - AMEND MANUALLY##\t")
                elif isinstance(elements[0], pmg.Element):
                    f.write(f"{elements[0].name}\t{elements[0].max_oxidation_state}\t{int(elements[0].name in CORE.LONE_PAIR_ELEMENTS)}\t")
                    osWarning = True
                elif isinstance(elements[0], pmg.Species):
                    f.write(f"{elements[0].element}\t{elements[0].oxi_state}\t{int(elements[0].element.name in CORE.LONE_PAIR_ELEMENTS)}\t")
                else:
                    raise Exception("Unexpected Site Contents")
            
            elif isinstance(site.species, pmg.Species):
                f.write(f"{site.species.symbol}\t{site.species.oxidation_state}\t{int(site.species.symbol in CORE.LONE_PAIR_ELEMENTS)}\t")
            else:
                f.write("##DISORDERED SITE - AMEND MANUALLY##\t")
            
            f.write(f"{site.coords[0]}\t{site.coords[1]}\t{site.coords[2]}\n")

        if osWarning:
            logging.warning("Oxidation States have been set to maximum (due to limitations with pymatgen). Amend input file with correct os")


    
# createInputFromCif("cif-files/1521543.cif", "files/na-abs.inp", "Na+")
# createInputFromCif("cif-files/Binary Fluorides/ICSD_CollCode5270 (beta-PbF2).cif", "files/na-abs.inp", "F-")
# bvDatToDb("cif-files/database_binary.dat", "soft-bv-params.sqlite3")
# ionDatToDb("cif-files/database_unitary.dat", "soft-bv-params.sqlite3")
# db = BVDatabase("soft-bv-params.sqlite3")
# print(db.getParams("Pb2+","O2-"))
        