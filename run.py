import logging, sys
from datetime import datetime
from bv2 import *
from fileIO import *

def create_input(args:list):

    if len(args) == 3:
        createInputFromCif(args[0], args[1], args[2])
    else:
        raise Exception(f"Method create_input requires 3 arguments, got {args}.\nThe following arguments are required: 'fileIn, 'fileOut', 'conductor'")

def bvs(args:list):

    if len(args) == 3:
        pbsnf4 = BVStructure.from_file(args[0])
        pbsnf4.initaliseMap(float(args[2]))
        pbsnf4.populateMapAbs()
        pbsnf4.exportMap(args[1], "delta")
    else:
        raise Exception(f"Method bvs requires 3 arguments, got {args}.\nThe following arguments are required: 'fileIn, 'fileOut', 'resolution'")
    
def bvs_penalty(args:list):

    if len(args) == 3 or len(args) == 5:
        pbsnf4 = BVStructure.from_file(args[0])
        pbsnf4.initaliseMap(float(args[2]))
        pbsnf4.createLonePairs()
        if len(args) == 3:
            pbsnf4.populateMismatchMap(penalty=0.05, fType="quadratic")
        else:
            pbsnf4.populateMismatchMap(penalty=float(args[3]), fType=args[4])
        pbsnf4.exportMap(args[1])
    else:
        raise Exception(f"Method bvs requires 3 arguments, got {args}.\nThe following arguments are required: 'fileIn, 'fileOut', 'resolution'")

def bvs_only_penalty(args:list):

    if len(args) == 3 or len(args) == 5:
        pbsnf4 = BVStructure.from_file(args[0])
        pbsnf4.initaliseMap(float(args[2]))
        pbsnf4.createLonePairs()
        if len(args) == 3:
            pbsnf4.populateMismatchMap(penalty=0.05, fType="quadratic", only_penalty=True)
        else:
            pbsnf4.populateMismatchMap(penalty=float(args[3]), fType=args[4], only_penalty=True)
        pbsnf4.exportMap(args[1])
    else:
        raise Exception(f"Method bvs requires 3 arguments, got {args}.\nThe following arguments are required: 'fileIn, 'fileOut', 'resolution'")

def site_bvs(args:list):
    if len(args) == 2:
        structure = BVStructure.from_file(args[0])
        structure.defineBufferArea()
        structure.findBufferedSites()
        
        for site in structure.sites.itertuples():
            print(f"Site {site.Index} at {site.coords} = {structure.findSiteBVS(site.Index, bool(int(args[1])))}")

def bulk_bvsmp(args:list):

    with open(args[0]) as f:
        for line in f.readlines():
            structure = BVStructure.from_file(line.strip())
            structure.initaliseMap(0.2)
            structure.createLonePairs()
            structure.dfToCif(f"{line.strip()[:-4]}-lp.cif")

            for k in [0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16]:
                for t in ["l", "q"]:
            # for k in [0.08]:
            #     for t in ["l"]:
                    structure.populateMismatchMap(penalty=k, fType=t)
                    structure.exportMap(f"{line.strip()[:-4]}-{k}{t}.grd")
                    structure.resetMap()

def render(args:list):

    if len(args) == 2 or len(args) == 3:
        structure = BVStructure.from_file(args[0])
        structure.defineBufferArea()
        structure.findBufferedSites()
        if len(args) == 3:
            if args[2] in ["lp", "lone-pair"]:
                structure.createLonePairs()
        logging.debug(structure.bufferedSites)
        structure.dfToCif(args[1])

def analyse(args:list):
    
    lessThan = np.vectorize(lambda x, y: float(x) < y, excluded=(1,))
    toInt = np.vectorize(lambda x: int(x))

    limits = [0.05, 0.075, 0.10]

    with open(args[0]) as f:
        contents = f.readlines()
        voxelNo = toInt(np.array(contents[2].split())).prod()
        
        print(f"Conductivity Volume Analysis for the file {args[0]}\nMismatch Level\tPercentage Under")
        for limit in limits:
            underLim = lessThan(np.array(contents[3].split()), limit)
            percentUnder = underLim.sum() / voxelNo * 100
            print(f"{limit}\t\t{percentUnder}")


# Code Allows Command Line Running of Functions https://stackoverflow.com/a/52837375
if __name__ == '__main__' and len(sys.argv) > 1:
    logging.basicConfig(level=logging.INFO, format='%(asctime)s -  %(levelname)s -  %(message)s', handlers=[
        # logging.FileHandler(f"./logs/{datetime.now().isoformat(timespec='seconds')}.log"),
        logging.FileHandler("logs/bvs.log"),
        logging.StreamHandler()
    ])

    start_time = datetime.now()

    if len(sys.argv) > 2:
        print(sys.argv)
        globals()[sys.argv[1]](sys.argv[2:])
    else:
        globals()[sys.argv[1]]()

    logging.info(f"Program Complete - Time Taken: {(datetime.now() - start_time)}")

else:

    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s -  %(levelname)s -  %(message)s', handlers=[
        # logging.FileHandler(f"./logs/{datetime.now().isoformat(timespec='seconds')}.log"),
        logging.StreamHandler()
    ])

    #analyse(["cif-files/bulk-bvsm/PbSnF4-0.08l.grd", 0.1])  
    # bvs_penalty(["files/pbsnf4.inp", "files/pbsnf4lp.grd", 0.25])
    render(["cif-files/bulk-bvsm/KSn2F5.inp","cif-files/bulk-bvsm/KSn2F5-gen.cif"])
    # create_input(["cif-files/ternary-fluorides/EntryWithCollCode152949 (PbSnF4).cif", "files/pbsnf4.inp", "F-"])




