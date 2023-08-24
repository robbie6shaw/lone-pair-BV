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
        pbsnf4.populateMap()
        pbsnf4.exportMap(args[1], "delta")
    else:
        raise Exception(f"Method bvs requires 3 arguments, got {args}.\nThe following arguments are required: 'fileIn, 'fileOut', 'resolution'")

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
    # bvs(["files/na-abs.inp", "files/na-abs-out.grd", 0.25])
    create_input(["cif-files/ternary-fluorides/EntryWithCollCode152949 (PbSnF4).cif", "files/pbsnf4.inp", "F-"])




