# dustfiles #

-Original Fortran files for Wiscombe's Mie scattering algorithms**, including test cases <br />
-Python translation of Wiscombe's ErrPack and MIEV0noP files <br />
-Draine's dust data tables <br />
-iPython notebook files plotting Mie and Draine values <br />
<br />

**Only edits made to the original Fortran files are print statements, to be used as checkpoints for debugging

Notes
-----
The Python translation of Wiscombe's MIEV0 file is included just in case, since a lot of work has already been done on it - however, it is definitely incomplete and needs a fair amount of editing/testing before it's usable.

Print checkpoints are not necessarily supposed to all show up or appear in the correct order; they're intended as tags that indicate where you are in the code when it's running.

Subroutines from the original Fortran code have been translated to functions in the Python code.
Definitions/specifications for all the Fortran subroutines were kept in the Python code in their entirety as comments.
