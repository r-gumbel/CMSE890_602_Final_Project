# **Benchmarking static calculations of nuclei using VU-TDHF3d**

##(_Pseudocode to outline general logic of workflow_)

### _Note: the emphasis here is on selecting pairs of nuclei which will later be used in a low energy (non-relativistic) collision.  This, then, provides the reason for selecting two nuclei at once, for their are certain decisions that need to prove consistent between pairs of nuclei (e.g. combining two ground state calculations with different interpolation setting will cause any time-dependent collision modeling to fail with this software)_

### 1.) User selects two nuclei of interest for they will use VU-TDHF3D to calculate their ground state properties.  Note: the emphasis here is on selecting pairs of nuclei which will later be used in a low energy (non-relativistic) collision. 
    * E.g. $$^{48}$$Ca (Z=20, N=28) and $$^{208}$$Pb (Z=82, N=126)