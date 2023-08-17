import pymatgen.core as pmg

def generate_structure(cifFile:str) -> pmg.Structure:
    return pmg.Structure.from_file(cifFile)

def generate_superstructure(structure:pmg.Structure) -> pmg.Structure:
    return structure.make_supercell((5,5,5), in_place=False)

def create_grid(structure:pmg.Structure, superstructure:pmg.Structure, a_res:int, b_res:int, c_res:int):

    # Find the vector to translate from one voxel to another
    lat = structure.lattice
    transA, transB, transC = (lat.a/a_res, lat.b/b_res, lat.c/c_res)

    #Define Grid List
    grid = []

    x,y,z = 0,0,0
    while x <= lat.a:
        while y <= lat.b:
            while z <= lat.c:
                grid.append((x,y,z))
                z = z + transC
            y = y + transB
        x = x + transA
