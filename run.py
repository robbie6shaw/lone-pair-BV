import logging
from datetime import datetime
from bv2 import *

logging.basicConfig(level=logging.INFO, format='%(asctime)s -  %(levelname)s -  %(message)s', handlers=[
    logging.FileHandler(f"./logs/{datetime.now().isoformat(timespec='seconds')}.log"),
    logging.StreamHandler()
    ])


pbsnf4 = BVStructure.from_file("pbsnf4.inp")
pbsnf4.initaliseMap(0.1)
pbsnf4.populateMap()
pbsnf4.exportMap("result.grd", "delta")