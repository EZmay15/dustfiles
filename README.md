# dustfiles #

-Original Fortran files for Wiscombe's Mie scattering algorithms**, including test cases 
-Python translation of the ErrPack and MIEV0noP files
-Draine's dust data tables
-iPython notebook files plotting Mie and Draine values

**Only edits made to the original files are print statements, to be used as checkpoints for debugging

Notes
-----
The MIEV0 file is included just in case, since a lot of it has already been translated to Python - however, it is definitely incomplete and needs a fair amount of editing/testing before it's usable.

Print checkpoints are not necessarily supposed to all show up or appear in the correct order; they are more like tags that indicate where you are in the code when it's running.

Subroutines from the original Fortran code have been translated to functions in the Python code.
Definitions/specifications for all the Fortran subroutines were kept in the Python code in their entirety as comments.
