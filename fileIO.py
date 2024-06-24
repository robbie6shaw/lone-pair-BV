import sqlite3, logging, re, collections
import CifFile as cf
import numpy as np
import pymatgen.core as pmg
from pathlib import Path          



class Ion:
    """
        Class that represents an ion as an object.
    """

    ION_REGEX = re.compile("([A-Za-z]{,2})(\d*)(\+|-)") # Regular expression for interpreting ions in string format
    LONE_PAIR_ELEMENTS = ["Pb", "Sn", "Bi", "Sb", "Tl"] # Elements that are searched for lone pairs
    one_reducer = lambda x: "1" if x == "" else x # Lambda function that removes the digit 1 from ions, e.g. Pb1+ -> Pb+

    def __init__(self, element:str, ox_state:int):
        """
            Initialises an ion class with the element and the oxidation state.
        """

        if isinstance(element, str):
            self.element = element
        else:
            raise TypeError("The element field of an Ion should be a string")
        
        self.ox_state = int(ox_state)
        self.string = self.element + (str(abs(self.ox_state)) if abs(self.ox_state) > 1 else "") + ("+" if self.ox_state > 0 else "-")
        self.radius = None

    def __str__(self):
        """
            Returns the string representation of the ion.
        """
        return self.string
    
    def __eq__(self, other):
        """
            Function to check whether the ion is equivalent to another ion, i.e. it is the same element and oxidation state
        """
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
        """
            Class method that allows creation of an ion object from the string representation.
            The string should be in the format of element then oxidation state, e.g. Pb2+, F-, Li+, O2-
        """
        result = cls.ION_REGEX.search(ion_string)
        
        if result is None:
            raise Exception(f"Cannot interpret the ion string sucessfully - {ion_string}")
        else:
            return Ion(element = result[1], ox_state = int(result[3] + cls.one_reducer(result[2])))
        
    def possible_lone_pair(self):
        """
            Checks whether the ion may have a lone pair
        """
        return self.element in self.LONE_PAIR_ELEMENTS



class BVDatabase:
    """
        A class representing a connection to a bond valence parameter database. Contains all methods required to communitate with the database.
    """
    
    DATABASE_DEFINITION = "database-define.sql"
    # Tuple for storing bond valence parameters
    bvparam = collections.namedtuple("BVParam", ['r0', 'ib', 'cn', 'r_cutoff', 'i1r', 'i2r', 'rmin', 'd0'])

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

        with open(self.DATABASE_DEFINITION, "r") as script:
            self.cursor.executescript(script.read())
        
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
            Creates a database entry for a BV parameter set. Only basic information, included in Brown's files is set in this method. Returns the row id.
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
        """
            Calculates expected value of the equilibrium bond distance R_min as defined by Chen et. al. 2019
        """
        x = (0.9185 + 0.2285 * abs(softness1 - softness2))*r0
        y = b*np.log(abs(osCation)/cn)
        return x - y
    
    def d0(self, b:float, os1:int, os2:int, block1:int, rmin:float, period1:int, period2:int):
        """
            Calculates the bond breaking energy D_0 as defined by Chen et. al. 2019
        """

        # Factor that depends on whether the cation is a s/p/d/f-block element
        if block1 <= 1: c = 1
        else: c = 2 

        return ((b**2)/2 * 14.4 * (c*(abs(os1*os2))**(1/c)) / (rmin*np.sqrt(period1 * period2)))
    
    def _params_error_check(self, ion1:Ion, ion2:Ion):
        """
            Fetches the result of get_bv_params and checks the result is well formed
        """
        result = self.fetch_all()
        if result is None or len(result) == 0:
            raise Exception(f"The combination of ions ({ion1}, {ion2}) are not on the BV Parameters Database")
        elif len(result) != 1:
            logging.warn(f"Multiple different database entries for the same ions ({ion1}, {ion2})")
        return result[0]

    def get_bv_params(self, ion1:Ion, ion2:Ion, bvse:bool = False):
        """
            Finds the parameters for a combination of two ions and returns them.

            Returns a bvparam named tuple, containing:
            
                'r0', 'ib', 'cn', 'r_cutoff', 'i1r', 'i2r', 'rmin', 'd0'
        """
        
        # If the ions are bonding
        if (ion1.ox_state * ion2.ox_state) < 0:

            # Retrieve all the information from the database
            #                    0   1  2   3   4         5         6         7            8            9          10         11        12
            self.execute("SELECT r0, b, ib, cn, r_cutoff, i1.radii, i2.radii, i1.softness, i2.softness, i1.period, i2.period, i1.block, i2.block FROM BVParam JOIN Ion i1 JOIN Ion i2 On BVParam.ion1 = i1.id AND BVParam.ion2 = i2.id WHERE (i1.symbol = ? AND i1.os = ? AND i2.symbol = ? AND i2.os = ?) OR (i2.symbol = ? AND i2.os = ? AND i1.symbol = ? AND i1.os = ?)", (ion1.element, ion1.ox_state, ion2.element, ion2.ox_state, ion1.element, ion1.ox_state, ion2.element, ion2.ox_state,))
            
            result = self._params_error_check(ion1, ion2)

            # If BVSE is not being done, the returned parameters are simpler
            if bvse:

                if ion1.ox_state > 0: 
                    cationOs = ion1.ox_state
                else: 
                    cationOs = ion2.ox_state
                    
                rmin = self.rmin(result[7], result[8], result[0], result[1], cationOs, result[3])
                d0 = self.d0(result[1], ion1.ox_state, ion2.ox_state, result[11], rmin, result[9], result[10])

                return self.bvparam(r0 = result[0], ib = result[2], cn = result[3], r_cutoff = result[4], i1r = result[5], i2r = result[6], rmin = rmin, d0 = d0)
            
            else:

                return self.bvparam(r0 = result[0], ib = result[2], cn = result[3], r_cutoff = result[4], i1r = None, i2r = None, rmin = None, d0 = None)

        # If BVSE and ions are repelling   
        elif bvse:

            return self.bvparam(r0=None, ib=None, cn=None, r_cutoff=None, i1r=self.get_radius(ion1), i2r=self.get_radius(ion2), rmin=None, d0=None)
        
        else:

            return None

        
    def get_atomic_no(self, element:str):
        """
            Function to find the atomic number of a particular element, given the symbol.
        """

        self.execute("SELECT atomic_no FROM Ion WHERE symbol = ? AND atomic_no IS NOT NULL", (element,))
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


def readCif(fileLocation:str|Path) -> cf.ReadCif:
    """
        Reads in a cif file from a Path object or a string of the path
    """
    return cf.ReadCif(str(fileLocation)).first_block()

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
    db.reset_database()

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

        started = False

        while True:

            entry = data.readline()
            
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
        

        started = False
        while True:
            
            entry = data.readline()

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


def create_input_from_cif(fileIn:Path, fileOut:Path, conductor:str):
    """
        Converts a crystal structure stored in a cif into a input file for further processing. Requires the input cif,
        the output file location and the conductor.
    """

    # Create a pymatgen object
    struct = pmg.Structure.from_file(fileIn)
    conductor = Ion.from_string(conductor)

    if fileIn.suffix != ".cif":
        logging.error(f"Incorrect file format. The input file should be a cif file and not a {fileIn.suffix} file")
    elif fileOut.suffix != ".inp":
        logging.error(f"Incorrect file format. The output file should be a inp file and not a {fileIn.suffix} file")

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
        disorderWarning = False

        # Create an entry in the dictionary to keep track of site multiplicity
        for site in struct.sites:
            siteDict[site.label] = 0

        # For every site, add label, element, os and cartesian coords
        for site in struct.sites:

            # Write the normal label and then a new label for P1 symmetry
            f.write(f"{site.label}\t{site.label}.{siteDict[site.label]}\t")
            siteDict[site.label] += 1

            if isinstance(site.species, pmg.Composition):
                
                elements = site.species.elements

                # If there are multiple elements on the site, throw warning.
                if len(elements) != 1:
                    f.write("##DISORDERED SITE - AMEND MANUALLY##\t")
                    disorderWarning = True
                
                # If only element for site defined, assume maximum oxidation state and throw warning.
                elif isinstance(elements[0], pmg.Element):
                    f.write(f"{elements[0].name}\t{elements[0].max_oxidation_state}\t{int(elements[0].name in Ion.LONE_PAIR_ELEMENTS)}\t")
                    osWarning = True

                # If species is defined for site, work normally
                elif isinstance(elements[0], pmg.Species):
                    f.write(f"{elements[0].element}\t{elements[0].oxi_state}\t{int(elements[0].element.name in Ion.LONE_PAIR_ELEMENTS)}\t")

                # Otherwise something very unexpected has happened    
                else:
                    raise Exception("Unexpected Site Contents")
            
            # If the site is only defined as a species, deal with that
            elif isinstance(site.species, pmg.Species):
                f.write(f"{site.species.symbol}\t{site.species.oxidation_state}\t{int(site.species.symbol in Ion.LONE_PAIR_ELEMENTS)}\t")

            else:
                f.write("##DISORDERED SITE - AMEND MANUALLY##\t")
                disorderWarning = True
            
            f.write(f"{site.coords[0]}\t{site.coords[1]}\t{site.coords[2]}\n")

        if osWarning:
            logging.warning("Oxidation States have been set to maximum (due to limitations with pymatgen). Amend input file with correct os")

        if disorderWarning:
            logging.warning(f"Strucutre has disordered sites - modification of ouput ({fileOut}) file is neeeded")



# TESTING AREA
    
# createInputFromCif("cif-files/1521543.cif", "files/na-abs.inp", "Na+")
# createInputFromCif("cif-files/Binary Fluorides/ICSD_CollCode5270 (beta-PbF2).cif", "files/na-abs.inp", "F-")
# bvDatToDb("cif-files/database_binary.dat", "soft-bv-params.sqlite3")
# ionDatToDb("cif-files/database_unitary.dat", "soft-bv-params.sqlite3")
# db = BVDatabase("soft-bv-params.sqlite3")
# print(db.getParams("Pb2+","O2-"))
        