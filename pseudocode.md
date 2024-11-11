# **Benchmarking static calculations of nuclei using VU-TDHF3d**

##(_Pseudocode to outline general logic of workflow_)

### _Note: the emphasis here is on selecting pairs of nuclei which will later be used in a low energy (non-relativistic) collision.  This, then, provides the reason for selecting two nuclei at once, for their are certain decisions that need to prove consistent between pairs of nuclei (e.g. combining two ground state calculations with different interpolation setting will cause any time-dependent collision modeling to fail with this software)_

### 1.) User selects two nuclei of interest for they will use VU-TDHF3D to calculate their ground state properties.  Note: the emphasis here is on selecting pairs of nuclei which will later be used in a low energy (non-relativistic) collision. 
    * E.g.<sup>48</sup>Ca (Z=20, N=28) and <sup>208</sup>Pb (Z=82, N=126)

### 2.) Prepare input files.  This consists mainly of the user selecting parameters in two main files--tdhf3d.inp and bspl.inp.
    * tdhf3d.inp: Describes entrance channel concerns--e.g. mass, charge, quadrupole/octupole contraints--as well as specifies key features of the computation itself.  It also specifies elements of the computation--e.g. convergence tolerance/max no. iterations for self-consistently solving differential equations.
    * bspl.inp: Sets dimensions of virtual box, as well as defines number of collocations/grid points for each each axes in our three dimensional cartesian coordinate space (x,y,z).

### 3.) Use makefile to build code.  Essentially takes this step to compile and link all of the various fortran procedures that comprise a given VU-TDHF3D run, which will run in parallel.  Makefile templates are available within the code for both licensed compilers (ifort) and open sourced (gfortran), and running any of these creates an executable _xtdhf_.

### 4.) Run job. Depending on such factors as computer hardware, total nuclear mass of the system to be calculated, dimensions of the virtural box, and certain user specified conditions for computation (e.g. max iterations and convergence tolerance), each static run can take anywhere from several minutes to hours.  

### 5.) While running, program writes to output file _tdhf3d.out_.  User can specify how often program writes to file with parameter _mprint_. (Standard setting is every 20 iterations of the program).  Program also generates a _.static_ file for each nucleus.

### 6.) These three files--_tdhf3d.out, A1.static, A2.static_--are of great interest to us.
    * _tdhf3d.out_ provides input to our benchmarking process.
    * _A1.static, A2.static_ provides output to be stored for later use in dynamic calculations ***provided*** the benchmarking in the previous step goes well.

### 7.) For each nuclei.
    * _EHF_: final ground state energy of system calculated with VU-TDHF3D.
    * _Beta_: Deformation of nucleus in its ground state.  This is calculated with the quadrupole moment associated with ground state energy.  

### 8.) Note: I'll need a section of the workflow devoted to extracting both of these g.s. observables. 
    *If both _EHF_ and _Beta_ are within an acceptable proximity to values calculated with Hartree-Fock-Bogoliubov; workflow is complete.
    * Else, jobs are flagged for one or both nuclei, error range is printed out, and user is prompted to return to step 2, restarting job with needed adjustments.  
    