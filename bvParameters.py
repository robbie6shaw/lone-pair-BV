import CifFile as cf
import sqlite3
import logging
import re

class BVDatabase:

    def __init__(self, dbLocation:str):
        """
            Initialise the connection to the database using its location.
        """
        self.conn = sqlite3.connect(dbLocation)
        self.cursor = self.conn.cursor()

    def commit(self):
        """
            Commit changes to the database
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
            Method for fetching all results from the database
        """
        return self.cursor.fetchall()

    def fetchone(self):
        """
            Method for fetching one result from the database
        """
        return self.cursor.fetchone()
    
    def extractFetchone(self):
        """
            Method for removing the tuple around the results from fetchone()
        """
        result = self.fetchone()
        if result is not None:
            return result[0]
        else:
            return None
    
    def lastRowId(self):
        """
            Method for getting the id of the last row changed on the database
        """
        return self.cursor.lastrowid

    def createDatabase(self):
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
                ib REAL
            );
            ''')
        
    def resetDatabase(self):
        """
            Method that resets the database
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
        self.execute("SELECT id FROM Ion WHERE symbol = ? AND os = ?", (ion, os, ))

        ans = self.extractFetchone()
        if ans == None:
            self.execute("INSERT OR IGNORE INTO Ion (symbol, os) VALUES (?,?)", (ion, os, ))
            return self.lastRowId()
        else:
            return ans

    def createEntry(self, ion1:str, os1:int, ion2:str, os2:int, r0:float, b:float):

        ion1Id = self.getOrInsertIon(ion1, os1)
        ion2Id = self.getOrInsertIon(ion2, os2)

        self.execute("SELECT id FROM BVParam WHERE ion1=? AND ion2=?", (ion1Id, ion2Id))
        if self.extractFetchone() == None:
            ib = 1/b
            self.execute("INSERT INTO BVParam (ion1, ion2, r0, b, ib) VALUES (?,?,?,?,?)", (ion1Id, ion2Id, r0, b, ib))

    def getParams(self, ion1, ion2):
        "F-""Pb2+""O2-"




def readCif(fileLocation:str):
    return cf.ReadCif(fileLocation).first_block()

def getDbCursor(dbLocation:str):
    conn = sqlite3.connect(dbLocation)
    return conn.cursor()

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

# db = BVDatabase("bv-params.sqlite3")
cifToDb("cif-files/bvparm2020.cif", "bv-params.sqlite3")
        