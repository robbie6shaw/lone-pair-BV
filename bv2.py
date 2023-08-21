from typing import Sequence
from numpy.typing import ArrayLike
import pymatgen.core as pmg
import math
import numpy as np
import pandas as pd
import fileIO

class BVStructure:

    DB_LOCATION = "soft-bv-params.sqlite3"

    def __init__(self, inputStr:str):
        """
            Initialises a BVStructure object using an input string.
        """
        
        lines = inputStr.splitlines()


        self.conductor = tuple(lines[0].split("\t"))
        self.params = tuple(lines[1].split("\t"))
        extraP = lines[2].split("\t")
        self.volume = extraP[0]
        self.vectors = np.zeros((3,3))
        for i in range(3,6):
            cols = lines[i].split("\t")
            for j in range(3):
                self.vectors[i-3][j] = cols[j]

        sites = []
        for i in range(7, len(lines)):
            data = lines[i].split("\t")
            sites.append({"label": data[0], "element": data[1], "ox_state": data[2], "a": data[3], "b": data[4], "c": data[5]})

        self.sites = pd.DataFrame(sites)

    def from_file(fileName:str):
        """
            Initialises a BVStructure object from an input file
        """
        with open(fileName, "r") as f:
            return BVStructure(f.read())
    
    def initaliseMap(self, resolution:int):
        """
            Initialises a map for storing the calculated BVS values. Creates a buffer cell structure, finds the core cells coordinates within that strcuture and defines the number of voxels. Arguments: \n
            resolution - Set a resolution for the map in armstrongs.
        """

        # Define the resolution and the intial values of the buffer area (ensures full coordination sphere is calculated)
        self.resolution = resolution
        self.bufferArea = np.array((3,3,3))



        # If one unit cell is less than the cutoff radius, include another unit cell to the buffer area
        for i in range(3):
            if self.coreCell.lattice.lengths[i] < self.rCutoff:
                self.bufferArea[i] += 2

        # Create a supercell from the buffer area using pymatgen
        self.bufferCell = self.coreCell.make_supercell(self.bufferArea, in_place=False)

        # Define the cartesian coordinates of the 'core cell' - the one that the BV calcualtions will take place in. This is achieved by find the cell number and then multiplying by the respective lattice parameter.
        self.coreCartesian = np.zeros(3)
        for i in range(3):
            self.coreCartesian[i] = math.floor(self.bufferArea[i]/2) * self.coreCell.lattice.lengths[i]

        # Calculate the number of voxels in each axis that is required to achieve the requested resolution
        self.voxelNumbers = np.zeros(3, dtype=int)

        for i in range(3):
            self.voxelNumbers[i] = math.ceil(self.coreCell.lattice.lengths[i] / self.resolution)

        # Initalise a map of dimensions that match the number of voxels
        self.map = np.zeros(self.voxelNumbers)

    def calcCartesian(self, shift:int):
        """
            Calculates the cartesian coordinates of voxel in the map using the origin of the 'core cell' and an integer shift. \n

            Returns a numpy array of floats defining the voxels cartesian coordinates.
        """
        position = np.copy(self.coreCartesian)
        
        for i in range(3):
            position[i] += shift[i] * self.resolution

        return position
    
    def calcDistanceWCutoff(self, point1, point2):
        """
            Calculates the distance between two points. If the the distance on one axis exceeds the cutoff distance, only that axes distance is returned.
        """

        deltaX = point2[0] - point1[0]
        if deltaX > self.rCutoff:  return deltaX
        deltaY = point2[1] - point1[1]
        if deltaY > self.rCutoff:  return deltaY
        deltaZ = point2[2] - point1[2]
        if deltaZ > self.rCutoff:  return deltaZ

        return math.sqrt(deltaX**2 + deltaY**2 + deltaZ**2)


    def populateMap(self, conductor:str):
        """
            Main function that populates the map space with BVS values. 
        """

        # Need to make this more general
        # Creates dataframe from bufferCell
        allAtoms = self.bufferCell.as_dataframe()

        print(self.coreCell)
        # Changes the composition row into a string representing the ion
        allAtoms["Species"] = allAtoms["Species"].map(lambda x: x.__str__())
        # Removes all conducting ions from the structure
        selectedAtoms = allAtoms[allAtoms["Species"] != conductor]

        print(self.coreCell)
        print(allAtoms["Species"])
        # Set up items to find bv parameters
        bvParams = {}
        db = fileIO.BVDatabase(self.DB_LOCATION)
        # For every ion that is not the conductor (currently assuming only one )
        for fixedIon in allAtoms["Species"]:
            if fixedIon == conductor:
                continue
            else:
                bvParams[fixedIon] = db.getParams(conductor, fixedIon)

        
        # For every voxel
        for h in range(self.voxelNumbers[0]):
            for k in range(self.voxelNumbers[1]):
                for l in range(self.voxelNumbers[2]):

                    # Calculate the voxels cartesian coordinates
                    pos = self.calcCartesian((h,k,l))

                    # Initialise the bond valence sum
                    bvSum = 0.

                    # For each atom in the structure
                    for i, fixedIon in selectedAtoms.iterrows():

                        # Calculate the point to point distance between the voxel position and the atom position
                        ri = self.calcDistanceWCutoff(pos, (fixedIon["x"], fixedIon["y"], fixedIon["z"]))

                        # If the seperation is less than 1 Å, set the BV value to very high value so the site is disregarded. This will cause the atom loop to be exited -> The site has a BV too high to be considered.
                        if ri < 1:
                            bvSum = 100
                            break

                        # If the seperation is greater than the cutoff radius, the bv contribution is 0.
                        elif ri > self.rCutoff:
                            continue

                        # Otherwise, calcualted the BV value and add it to the total
                        else:
                            r0, ib = self.BV_PARAMS[fixedIon["Species"]]
                            bv = calcBV(r0, ri, ib)
                            bvSum += bv

                    # Update the map
                    self.map[h][k][l] = bvSum

    def deltaBV(self, value:float, ion:str):
        if ion == "F-" or ion == "Na+":
            result = abs(value - 1)
            return result

    def exportMap(self, fileName:str, dataType:str):

        if dataType == "delta":
            deltaBV_vector = np.vectorize(self.deltaBV, excluded=("ion"))
            export = deltaBV_vector(self.map, "F-")
        else:
            export = self.map

        print(export)

        with open(fileName, 'w') as file:
            file.write("%s\n" % (self.coreCell.formula))
            file.write("%f %f %f %f %f %f\n" % self.coreCell.lattice.parameters)
            file.write("%i %i %i\n" % tuple(self.voxelNumbers.tolist()))

            export.tofile(file,"  ")


def generate_structure(cifFile:str) -> pmg.Structure:
    return pmg.Structure.from_file(cifFile)

def calcBV(r0:float, ri:float, ib:float) -> float :
    """
        Calculate the bond valence from distance. Arguments: \n
        r0 - The radius bond valence parameter \n
        ri - The current distance \n
        ib - The inverse of the bond valence parameter \n

        Will be implemented in fortran to improve performance
    """
    return math.exp((r0 - ri) * ib)

def calcPPDistance(point1, point2):
    deltaX = point2[0] - point1[0]
    deltaY = point2[1] - point1[1]
    deltaZ = point2[2] - point1[2]

    return math.sqrt(deltaX**2 + deltaY**2 + deltaZ**2)

def calcPSDistance(x:float, y:float, z:float, site:pmg.PeriodicSite) -> float:
    """
        Calculate the distance between a point and a site. Arguments: \n
        x,y,z - Cartesian coordinates of the point \n
        site - A Pymatgen Object representing the site
    """
    return site.distance_from_point((x,y,z))

def calcSSDistance(site1:pmg.PeriodicSite, site2:pmg.PeriodicSite) -> float:
    """
        Calculates the distance between two sites. Does not use the pmg.PeriodicSite.distance as it does strange things with periodicity.
    """
    return calcPSDistance(site2.x, site2.y, site2.z, site1)

def findSiteBVS(site:pmg.PeriodicSite, structure:pmg.Structure, radius:int = 6 ) -> float:
    """
        Calculates the bond valence sum for a particular site located in a structure. The stucture must be ordered. Arguments: \n
        site - The site to calculate the bvs from. \n
        structure - The structure that the site is located in. \n
        radius - The cutoff radius for the bvs calculation. Defaults to 6 Å. \n
    """

    bvParams = {'Sn2+':(1.925, 2.702), 'Pb2+':(2.03, 2.702)}
    
    # Create temporary copy of structure
    tempStruct = structure.copy()

    # Ensure the site chosen only has one element on it
    if len(site.species.elements) != 1:
        raise Exception("Site chosen for BVS sum is disordered - must choose site that is ordered")

    # Find the species within the defined radius
    coordAtoms = structure.get_neighbors(site, radius)

    # Set bond valence sum to 0
    bvs = 0.

    for atom in coordAtoms:
        
        # If the oxidation state is the either both positive or both negative, disregard it from the bond valence sum
        # bool1 = atom.specie.oxi_state > 0
        # bool2 = site.specie.oxi_state > 0
        # bool3 = bool1 ^ bool2
        if ((atom.specie.oxi_state > 0) ^ (site.specie.oxi_state > 0)):

            # print(atom, end=' : ')
            
            # Find the bond valence parameters
            r0, ib = bvParams[atom.species_string]

            # Add this contribution to the bond valence
            bv = calcBV(r0, calcSSDistance(site, atom), ib)
            # print(calcSSDistance(site, atom), end=", ")
            # print(bv)
            bvs += bv
        
    return bvs

def bvsCif (fileLocation:str):
    """
        Find bond valence sum for all fluoride sites in a structure.
    """
    structure = generate_structure(fileLocation)

    for site in structure.sites:
        if site.species_string == "F-":
            print(f"F- Site at ({site.x}, {site.y}, {site.z}): {findSiteBVS(site, structure)}")

# bvsCif("cif-files/Binary Fluorides/ICSD_CollCode5270 (beta-PbF2).cif")
pbsnf4 = BVStructure.from_file("pbsnf4.inp")
print(pbsnf4.sites)
# pbsnf4.initaliseMap(1)
# pbsnf4.populateMap("Na+")
# pbsnf4.exportMap("result.grd", "delta")