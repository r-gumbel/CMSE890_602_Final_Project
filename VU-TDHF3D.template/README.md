# VU-TDHF3D
### INSTRUCTIONS FOR VU-TDHF3D

Input is read from three files:

- bspl.inp   -  defines the mesh
- tdhf3d.inp -  all other input variables
- skyrme.inp -  different skyrme forces
- pair.inp   -  pairing parameters

See top of subroutine "getin.f90" for all the input options and
files.

Main subroutine is tdhf3d.f90

Output is in file tdhf3d.out and tdhf3d.plot (plotmtv file).

The code is compiled using the shell scripts "build\_\*" in the top level
directory of the package. There is a different build for different
compilers, calling different makefiles in directories "bslib90" and
"tdhf".

To use a new different compiler, copy some of the makefiles to another
name and adjust the flags in the files for the appropriate compiler.
The makefiles with "omp" link with OpenMP libraries, whereas the
serial make files use the appropriate dummy library.
Before recompiling issue the cleaning option
```Bash
clean
```
which removes all the objects etc. from the previous compilation.

RUNNING:

The executable "xtdhf" is copied to the "run" subdirectory by
default. To run the executable using nohup do;
```Bash
nohup run &
```
the run script also calls the time command.

For any questions contact us.
