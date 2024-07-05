import logging, sys
from argparse import ArgumentParser
from datetime import datetime
from bvStructure import *
from fileIO import *
from pathlib import Path
from shutil import copy2

RESOLUTION_ARGS = {'default':0.1, 'type':float, 'help':"The target resolution of the produced map. The number of voxels will be rounded up to ensure divsibility by 12. Defaults to 0.1."}
EC_ARGS = {'action':'store_true', 'help':'Toggles whether the effective or absolute charge is used for repulsion calculations in bond valence site energy. Defaults to absolute charge.'}
NJ_ARGS = {'action':'store_true', 'help':'Toggles whether just-in-time compliation is used in the calcualtion. Defaults to using JIT for large speed gains, flag turns it off.'}

def create_input(parser:ArgumentParser, overrideArgs:list = None):

    parser.add_argument("cif_file", help="The structure to be analysed, in a cif file format.")
    parser.add_argument("output_location", help="The location for the create input file to be outputted to.")
    parser.add_argument("conductor", help="The conducting ion under investigation. Specified in the format (ELEMENT)(CHARGE NUMBER)(CHARGE SIGN)")
    args = parser.parse_args(overrideArgs)

    create_input_from_cif(Path(args.cif_file), Path(args.output_location), args.conductor)

def bvsm(parser:ArgumentParser, overrideArgs:list = None):

    parser.add_argument("input_file")
    parser.add_argument("output_file")
    parser.add_argument("-r", "--resolution", **RESOLUTION_ARGS)
    parser.add_argument("-m", "--mode", default=1, choices=range(0,3), type=int)
    parser.add_argument("-n", "--no_jit", **NJ_ARGS)
    parser.add_argument("-k", "--penalty_constant", default=0.05, type=float)
    parser.add_argument("-t", "--penalty_type", default="q", choices=("q","l","quadratic","linear"))

    args = vars(parser.parse_args(overrideArgs))
    args.pop('function', None)
    _bvsm(**args)

def _bvsm(input_file:str, output_file:str, resolution:float, mode:int, no_jit:bool, penalty_constant:float, penalty_type:str):

    crystal = BVStructure.from_file(input_file, bvse=True)
    crystal.initalise_map(resolution)

    if mode > 0:
        crystal.create_lone_pairs()

    if no_jit:
        if mode == 0:
            crystal.populate_map_bvsm(fType=penalty_type, penalty=0)
        elif mode == 1:
            crystal.populate_map_bvsm(fType=penalty_type, penalty=penalty_constant)
        elif mode == 2:
            crystal.populate_map_bvsm(fType=penalty_type, penalty=penalty_constant, only_penalty=True)

    else:
        if penalty_type == "l" or penalty_type == "linear":
            logging.error("Linear penalty functions are not implemented using JIT. Add flag --no_jit to run.")

        crystal.populate_map_bvsm_jit(mode = mode, penalty=penalty_constant)

    crystal.export_map(output_file)

def bvse(parser:ArgumentParser, overrideArgs:list = None):

    parser.add_argument("input_file")
    parser.add_argument("output_file")
    parser.add_argument("-r", "--resolution", **RESOLUTION_ARGS)
    parser.add_argument("-m", "--mode", default=1, choices=range(0,3), type=int)
    parser.add_argument("-e", "--effective_charge", **EC_ARGS)
    parser.add_argument("-n", "--no_jit", **NJ_ARGS)

    args = vars(parser.parse_args(overrideArgs))
    args.pop('function', None)
    _bvse(**args)

def _bvse(input_file:str, output_file:str, resolution:float, mode:int, effective_charge:bool, no_jit:bool):

    crystal = BVStructure.from_file(input_file, bvse=True)
    crystal.initalise_map(resolution)
    if mode > 0:
        crystal.create_lone_pairs()
    if no_jit:
        crystal.populate_map_bvse(mode=mode)
    else:
        crystal.populate_map_bvse_jit(mode = mode, effectiveCharge=effective_charge)
    crystal.export_map(output_file)


def site_bvs(parser:ArgumentParser, overrideArgs:list = None):

    parser.add_argument("input_file")
    parser.add_argument("-v", "--vector", action="store_true")

    args = parser.parse_args(overrideArgs)

    crystal = BVStructure.from_file(args.input_file)
    crystal.define_buffer_area()
    crystal.find_buffer_sites()
        
    for site in crystal.sites.itertuples():
        print(f"Site {site.Index} at {site.coords} = {crystal.find_site_bvs(site.Index, args.vector)}")

def _create_dir(path:Path, name:str):
    resultPath = path.joinpath(name)
    if not resultPath.is_dir():
        resultPath.mkdir()
    return resultPath

def bulk_bvse(parser:ArgumentParser, overrideArgs:list = None):

    parser.add_argument("base_path", help="The folder that the program should use for the calculations. Should contain a folder named 'cif' that contains all the structures to process. The results will be outputted to 'result.")
    parser.add_argument("conductor", help="The conducting ion under investigation. Specified in the format (ELEMENT)(CHARGE NUMBER)(CHARGE SIGN)")
    parser.add_argument("-r", "--resolution", **RESOLUTION_ARGS)
    parser.add_argument("-e", "--effective_charge", **EC_ARGS)
    args = parser.parse_args(overrideArgs)

    basePath = Path(args.base_path)
    cifPath = basePath.joinpath("cif")

    if not cifPath.is_dir():
        logging.error("The folder specified does not contain a folder called 'cif'")
        sys.exit()

    else:

        resultPath = basePath.joinpath("result")
        if resultPath.is_dir():
            for i in range(100):
                trialPath = resultPath.with_stem(f"{resultPath.stem}-{i}")
                if not trialPath.is_dir():
                    resultPath = trialPath
                    break

        resultPath.mkdir()

        for cifFile in cifPath.iterdir():
            formula = readCif(cifFile)['_chemical_formula_structural'].replace(' ','')
            # create_input_from_cif(cifFile, inpPath.joinpath(cifFile.name).with_suffix(".inp"), conductor)
            formulaFolder = _create_dir(resultPath, formula)
            inpFile = formulaFolder.joinpath(formula).with_suffix(".inp")
            cubeFile = formulaFolder.joinpath(formula).with_suffix(".cube")

            try:

                copy2(cifFile, formulaFolder.joinpath(formula).with_suffix(".cif"))

                create_input_from_cif(cifFile, inpFile, args.conductor)

                _bvse(input_file=inpFile, output_file=cubeFile, resolution=args.resolution, mode=1, effective_charge=args.effective_charge, no_jit=False)
            
            except Exception as e:
                logging.error(f"The following {type(e)} exception was raised when processing the structure {formula}. The following traceback was produced:")
                print(e)
                continue

def render(parser:ArgumentParser, overrideArgs:list = None):

    parser.add_argument("input_file")
    parser.add_argument("output_file")
    parser.add_argument("-l","--lp","--lone_pair", action='store_true', help='Toggles whether lone pairs are rendered in the cif file by adding dummy helium sites')
    args = parser.parse_args(overrideArgs)

    structure = BVStructure.from_file(args.input_file)
    structure.define_buffer_area()
    structure.find_buffer_sites()
    if args.lp:
        structure.create_lone_pairs()
    logging.debug(structure.bufferedSites)
    structure.export_cif(args.output_file)


def data_import(parser:ArgumentParser, overrideArgs:list = None):

    bvDatToDb("cif-files/database_binary.dat", "soft-bv-params.sqlite3")
    ionDatToDb("cif-files/database_unitary.dat", "soft-bv-params.sqlite3")

def buffer_export(parser:ArgumentParser, overrideArgs:list = None):
    """
        Export the buffered site list to an excel file.
    """
    
    parser.add_argument("input_file")
    parser.add_argument("output_file")
    args = parser.parse_args(overrideArgs)

    pbsnf4 = BVStructure.from_file(args.input_file)
    pbsnf4.initalise_map(1.0)
    pbsnf4.bufferedSites.to_excel(args.output_file)

# Code Allows Command Line Running of Functions https://stackoverflow.com/a/52837375
if __name__ == '__main__' and len(sys.argv) > 1:
    logging.basicConfig(level=logging.INFO, format='%(asctime)s -  %(levelname)s -  %(message)s', handlers=[
        # logging.FileHandler(f"./logs/{datetime.now().isoformat(timespec='seconds')}.log"),
        logging.FileHandler("logs/bvs.log"),
        logging.StreamHandler()
    ])

    start_time = datetime.now()

    parser = ArgumentParser(
        prog="lpbv",
        description="Tool that can calculate bond valence mismatch or bond valence site energy maps of a crystal structure using a cif as input. Can add dummy lone pair sites to the structure."
    )
    parser.add_argument(
        'function', 
        help="Specifies the function of the lpbv to use. Consult documentation for list of possible calls"
    )

    print(sys.argv)
    try:
        globals()[sys.argv[1]](parser)
        logging.info(f"Program Complete - Time Taken: {(datetime.now() - start_time)}")
    except KeyError:
        print("Invalid Function Entered. Possible options: create_input, bvsm, bvse, site_bvs, bulk_bvse, render, data_import, buffer_export")

else:

    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s -  %(levelname)s -  %(message)s', handlers=[
        # logging.FileHandler(f"./logs/{datetime.now().isoformat(timespec='seconds')}.log"),
        logging.StreamHandler()
    ])

    parser = ArgumentParser()

    # bvs(["files/pbsnf4.inp", "files/temp.grd", 1])
    # bvse_jit(['results/PbSnF4/PbSnF4.inp', 'results/PbSnF4/PbSnF4-test.cube', '1', '1'])
    # render(['results/Sn-II-database/result/Sn3BrF5/Sn3BrF5.inp', 'results/Sn-II-database/result/Sn3BrF5/Sn3BrF5-lp.cif','lp'])
    # create_input(["cif-files/ternary-fluorides/EntryWithCollCode152949 (PbSnF4).cif", "files/pbsnf4.inp", "F-"])
    # new_bulk(["/home/rs/bv-project/results/Sn-II-database", "F-", 0.1])
    # bvse(parser, ['results/Pb-II-database/result/PbPdF4/PbPdF4.inp', 'results/Pb-II-database/result/PbPdF4/PbPdF4.cube'])
    #bvse(parser, ['results/PbSnF4/PbSnF4.inp', 'results/PbSnF4/refactor-test/bvse06.cube'])
    #site_bvs(parser, ['results/PbF2-beta/PbF2.inp'])