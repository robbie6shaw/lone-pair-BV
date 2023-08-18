import pymatgen.core as pmg
import math
import numpy as np

def generate_structure(cifFile:str) -> pmg.Structure:
    return pmg.Structure.from_file(cifFile)

def create_centred_grid(structure:pmg.Structure, enlargementSize:int, resolution:int):


    enlargedcell = structure.make_supercell((enlargementSize, enlargementSize, enlargementSize), in_place=False)

    # Find the vector to translate from one voxel to another
    lat = structure.lattice
    transA, transB, transC = (lat.a/resolution, lat.b/resolution, lat.c/resolution)


    # Find the 0,0,0 coordinate of the examined cell
    middle = int(enlargementSize/2) + 1
    startCoords = (lat.a*middle, lat.b*middle, lat.c*middle)

    #Define Grid List
    grid = ()

    #Iterate through 3D space to find all points for the grid
    x,y,z = startCoords
    while x <= startCoords[0] + lat.a:

        while y <= startCoords[1] + lat.b:

            while z <= startCoords[2] + lat.c:
                grid.append((x,y,z))
                z = z + transC

            y = y + transB
            z = startCoords[2]

        x = x + transA
        y = startCoords[1]

    return enlargedcell, grid

def find_neigbours(eCell:pmg.Structure, grid:[(int,int,int)], radius:float):
    result = {}

    for pixel in grid:
        result[pixel] = eCell.get_sites_in_sphere(pixel, radius)

    return result

# def find_distance(site1, site2):
#     deltaX = site2[0] - site1[0]
#     deltaY = site2[1] - site1[1]
#     deltaZ = site2[2] - site1[2]

#     return math.sqrt(deltaX^2 + deltaY^2 + deltaZ^2)

def find_distances_to_neigbours(eCell:pmg.Structure, grid:[(int,int,int)], radius:float):
    result = {}

    for pixel in grid:
        neighbours = eCell.get_sites_in_sphere(pixel, radius)
        pixelList = []
        for neighbour in neighbours:
            pixelList.append((neighbour.species_string, neighbour.distance_from_point(pixel)))
        result[pixel] = pixelList

    return result

def bv_equation(r0, ri, b):
    return math.exp((r0-ri)/b)

def calculate_bv(distanceList:dict, bvParams):

    result = {}
    
    for pixel, data in distanceList.items():
        
        bvSum = 0

        for site, distance in data:
            bvSum += bv_equation(bvParams[site][0], distance, bvParams[site][1])

        result[pixel] = bvSum

    return result
        

def export_grd(fileName, result, res, structure:pmg.Structure):
    with open(fileName, 'w') as file:
        file.write("%s\n" % (structure.formula))
        file.write("%f %f %f %f %f %f\n" % structure.lattice.parameters)
        file.write("%i %i %i\n" % (res, res, res))

        lnProgress = 0
        for value in result.values():
            file.write("    %.10e" % (value))
            lnProgress += 1
            if lnProgress == 6:
                file.write("\n")
                lnProgress = 0


# resolution = 20
# bvParams = {'Sn2+':(1.925, 0.37), 'Pb2+':(2.03, 0.37)}
# structure = generate_structure("cif-files/Ternary Fluorides/EntryWithCollCode152949 (PbSnF4).cif")
# structure.remove_species(["F-"])
# print(structure)
# eCell, grid = create_centred_grid(structure, 5, resolution)
# distanceList = find_distances_to_neigbours(eCell, grid, 6)
# result = calculate_bv(distanceList, bvParams)
# export_grd("PbSnF2.grd", result, resolution, structure)
