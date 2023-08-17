import pymatgen.core as pmg

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
    grid = []

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


structure = generate_structure("cif-files/Ternary Fluorides/EntryWithCollCode152949 (PbSnF4).cif")
eCell, grid = create_centred_grid(structure, 5, 5)
print(find_neigbours(eCell, grid, 5))