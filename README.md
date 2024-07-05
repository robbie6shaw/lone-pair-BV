# lone-pair-BV

## Running The Software

To use lone-pair-BV, download the repository. The dependencies then have to be installed for which use of mambaforge is recommended for fast installation. To install the dependencies, the following command should be entered in the install folder:

```
  mamba env create -f .\conda-env.yml
```

Once the environment is activated, lone-pair-BV is ready to use. Use of the program consists of running the python script on the command line. Once in the install directory, `run.py` can run the following commands:

### create_input
The create input command is used to take a cif file and create a human readable + editable input file for creating a bond valence map. Positional arguments:
-  `cif_file` -        The structure to be analysed, in a cif file format.
-  `output_location` -  The location for the create input file to be outputted to.
-  `conductor` -        The conducting ion under investigation. Specified in the format
                   (ELEMENT)(CHARGE NUMBER)(CHARGE SIGN)
   
Options:
-   `-h, --help` -       show this help message and exit

### bvsm
The bond valence sum mismatch command creates a bond valence mismatch map for the structure. It accepts three arguments:
- The location of the input file
- The desired location of the output file
- The resolution of the voxels in armstrongs

### bvs_penalty
The bond valence sum with penalty command creates a bond valence mismatch map for the structure. Unlike the `bvs` command, it does apply a penalty function and creates dummy lone pair sites on some heavy metal atoms. 
It accepts the same arguments as the `bvs` command.

### site_bvs
A command to find the bond valence parameters for every site in a crystal structure. It accepts two arguments:
- The location of the input file
- An integer representing whether a vector or scalar sum is required. Enter 0 for scalar; enter 1 for vector.
