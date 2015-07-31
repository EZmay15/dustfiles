

#Tests MIEV0noP for debugging - may need some editing itself!

from Wis_MIEV0noP import *
import numpy as np


#PROGRAM TSTMIV

"""
c         Test cases for MIEV0 from NCAR Tech Note, with, in addition,
c         calculations of Legendre coefficients

c     NOTE:  set NoPMOM = True if using NoPMOM version of MIEV0


      IMPLICIT  NONE

c ----------------------------------------------------------------------
c --------  I / O SPECIFICATIONS FOR SUBROUTINE MIEV0  -----------------
c ----------------------------------------------------------------------
      INTEGER  MAXANG, MOMDIM
      PARAMETER  ( MAXANG = 40, MOMDIM = 10 )
      LOGICAL  ANYANG, PERFCT, PRNT(2)
      INTEGER  IPOLZN( 4,2 ), NUMANG, NMOM
      REAL     GQSC, MIMCUT, PMOM( 0:MOMDIM, 4 ), QEXT, QSCA, SPIKE,
     $         XMU(MAXANG), XX(4)
      COMPLEX  CREFIN( 2 ), SFORW, SBACK, S1( MAXANG ), S2( MAXANG ),
     $         TFORW( 2 ), TBACK( 2 )
c ----------------------------------------------------------------------

c     .. Local Scalars ..

      LOGICAL   NOPMOM
      INTEGER   I, NIOR, NXX
      REAL      PI
c     ..
c     .. External Subroutines ..

      EXTERNAL  MIEV0
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ASIN, COS
c     ..
      DATA  CREFIN / (1.5,0.0), (1.5,-0.1) /,
     $      XX / 10., 100., 1000., 5000. /,
     $      IPOLZN / 13, 24, 1234, 0, -13, -24, -1234, 0 /,
     $      PRNT / 2*.TRUE. /,
     $      PERFCT /.FALSE./,
     $      MIMCUT / 1.E-6 /,
     $      ANYANG /.TRUE./,
     $      NUMANG / 19 /
"""

MAXANG = 40
MOMDIM = 10
XMU = []
for i in range(MAXANG):
    XMU.append(0)
    
S1 = []
for i in range(MAXANG):
    S1.append(0)

S2 = []
for i in range(MAXANG):
    S2.append(0)
    
TFORW = [0,0]
TBACK = [0,0]

CREFIN = [ 1.5 + 0.0j, 1.5 - 0.1j ]
XX = [ 10., 100., 1000., 5000. ]
IPOLZN = np.zeros((4,2))
IPOLZN[0,0] = 13
IPOLZN[1,0] = 24
IPOLZN[2,0] = 1234
IPOLZN[3,0] = 0
IPOLZN[0,1] = -13
IPOLZN[1,1] = -24
IPOLZN[2,1] = -1234
IPOLZN[3,1] = 0
PRNT = [ True, True ]
PERFCT = False
MIMCUT = 1.E-6
ANYANG = True
NUMANG = 19
QEXT = None
QSCA = None
GQSC = None
PMOM = [1,2,3,4,5,6,7,8,9,10]
SFORW = None
SBACK = None
SPIKE = None

NOPMOM = True
PI     = 2.*np.arcsin( 1.0 )

if NUMANG > MAXANG:
    print "numang"
    import sys
    sys.exit()

for I in range( 0, NUMANG ):
    XMU[ I ] = np.cos( (I)*np.pi / ( NUMANG - 1 ) )

XMU[ NUMANG - 1 ] = -1.0


for NIOR in range( 0, 2 ):

    for NXX in range( 0, 4 ):

        NMOM = 10
        if NOPMOM or NXX == 3:
            NMOM = 0

        print MIEV0( XX[ NXX ], CREFIN[ NIOR ], PERFCT, MIMCUT, ANYANG, NUMANG,
                    XMU, NMOM, IPOLZN[ NXX, NIOR ], MOMDIM, PRNT, QEXT, QSCA,
                    GQSC, PMOM, SFORW, SBACK, S1, S2, TFORW, TBACK, SPIKE )

