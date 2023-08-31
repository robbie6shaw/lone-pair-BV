# bv-project

## Running The Software

To utilise bv-project, it is run using python on the command line. Once in the 'bv-project' directory, `run.py` can run the following commands:

### create_input
The create input command is used to take a cif file and create a human readable + editable input file for creating a bond valence map. It accepts three arguments:
- The location of the cif file
- The desired location of the output file
- The ion that will conduct

## bvs 
The bond valence sum command creates a bond valence mismatch map for the structure. It does not apply any penalty functions or create lone pairs. It accepts three arguments:
- The location of the input file
- The desired location of the output file
- The resolution of the voxels in armstrongs

## bvs_penalty
The bond valence sum with penalty command creates a bond valence mismatch map for the structure. Unlike the `bvs` command, it does apply a penalty function and creates dummy lone pair sites on some heavy metal atoms. 
It accepts the same arguments as the `bvs` command.

## site_bvs
A command to find the bond valence parameters for every site in a crystal structure. It accepts two arguments:
- The location of the input file
- An integer representing whether a vector or scalar sum is required. Enter 0 for scalar; enter 1 for vector.
