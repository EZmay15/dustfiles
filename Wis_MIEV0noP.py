

#The original Wiscombe MIEV0noP Fortran file is also in the dustfiles repository for
#reference/comparison as MIEV0noPcopy.f! MIEV.doc, ErrPack.f, and MIEV0.f have also been
#included as copies. I found the code on the page "Mie Scattering" by Scott Prahl, who had
#zipped a version of Wiscombe's code and made it downloadable (http://omlc.org/software/mie/).

#Notes:
#Most for loops in this code loop through calling elements of arrays. Because Fortran
#arrays are indexed starting from 1 while Python arrays are indexed from 0, I changed
#many (but not all) loops to begin from 0 instead of 1. Many references to list elements are
#also corrected in this way (ex. XMU[1] -> XMU[0]). Operations that would be affected by these
#changes have been edited as well (often by adding or subtracting 1).
#Complex numbers have been changed from Fortran's (0., 0.) format to Python's 0. + 0.j format
#Definitions/specifications for subroutines in the original Fortran file
#(subroutines -> functions in this Python code) were kept in this code as comments.
#If output is formatted oddly, the first place to look in this code would be the print
#statements, since those have been converted in the simplest way possible and may have lost
#some specifications possible with Fortran's "write" statement.


import numpy as np


def MIEV0( XX, CREFIN, PERFCT, MIMCUT, ANYANG, NUMANG, XMU, NMOM, IPOLZN, MOMDIM, 
           PRNT, QEXT, QSCA, GQSC, PMOM, SFORW, SBACK, S1, S2, TFORW, TBACK, SPIKE ):
    
    """
        NoPMOM version:  saves storage by deleting everything
        having to do with PMOM calculation:

        * deletes modules LPCoef, LPCO1T, LPCO2T

        * NMOM,IPOLZN,PMOM,MOMDIM,CalcMo,NPQuan deleted from argument
          lists and everywhere else they occur internally, although the
          first four are retained in the MIEV0 argument list for
          compatibility with complete version of MIEV0

        * MIEV0:   shrinks size of LitA,LitB arrays to 2;  deletes call
                   to LPCoef;  saving of An,Bn in LitA,LitB arrays deleted

        * MiPrnt:  PMOM print deleted

        * CkInMi:  NMOM test changed to allow only zero value (nonzero
                   value would have no impact but indicates user
                   conceptual error or negligence)

        * TestMi:  saving and restoring of NMOM,IPOLZN deleted;
                   30-loop deleted
                   
    Computes Mie scattering and extinction efficiencies; asymmetry
    factor;  forward- and backscatter amplitude;  scattering
    amplitudes vs. scattering angle for incident polarization parallel
    and perpendicular to the plane of scattering;
    some quantities needed in polarized radiative transfer;  and
    information about whether or not a resonance has been hit.
    
    Input and output variables are described in file MIEV.DOC;
    internal variables are described in MIEV0.f file.
    
    CALLING TREE:
    
        MIEV0
            TESTMI
                TSTBAD
                MIPRNT
                ERRMSG
            CKINMI
                WRTBAD
                WRTDIM
                ERRMSG
            SMALL1
            SMALL2
            ERRMSG
            BIGA
                CONFRA
                    ERRMSG
            MIPRNT


        IMPLICIT NONE
    
    ----------------------------------------------------------------------
    --------  I / O SPECIFICATIONS FOR SUBROUTINE MIEV0  -----------------
    ----------------------------------------------------------------------
        LOGICAL  ANYANG, PERFCT, PRNT(*)
        INTEGER  IPOLZN, MOMDIM, NUMANG, NMOM
        REAL     GQSC, MIMCUT, PMOM( 0:MOMDIM, * ), QEXT, QSCA, SPIKE, XMU(*), XX
        COMPLEX  CREFIN, SFORW, SBACK, S1(*), S2(*), TFORW(*), TBACK(*)
    ----------------------------------------------------------------------
    
    ** NOTE -- MaxTrm = 10100 is necessary to do some of the test probs,
               but 1100 is sufficient for most conceivable applications
    
    .. Parameters ..
    
    INTEGER   MAXANG, MXANG2
    PARAMETER ( MAXANG = 501, MXANG2 = MAXANG / 2 + 1 )
    INTEGER   MAXTRM
    PARAMETER ( MAXTRM = 10100 )
    REAL      ONETHR
    PARAMETER ( ONETHR = 1. / 3. )
    ..
    .. Local Scalars ..
    
    LOGICAL   NOABS, PASS1, YESANG
    INTEGER   I, J, N, NANGD2, NTRM
    REAL      CHIN, CHINM1, COEFF, DENAN, DENBN, FN, MIM, MM, MRE, NP1DN, PSIN, 
              PSINM1, RATIO, RIORIV, RN, RTMP, TAUN, TCOEF, TWONP1, XINV
    COMPLEX   AN, ANM1, ANP, ANPM, BN, BNM1, BNP, BNPM, CDENAN, CDENBN, CIOR, 
              CIORIV, CTMP, ZET, ZETN, ZETNM1
    ..
    .. Local Arrays ..
        
    REAL      PIN( MAXANG ), PINM1( MAXANG ), RBIGA( MAXTRM )
    COMPLEX   CBIGA( MAXTRM ), LITA( 2 ), LITB( 2 ), 
              SM( MAXANG ), SMS( MXANG2 ), SP( MAXANG ), SPS( MXANG2 )
    ..
    .. External Subroutines ..
    
    EXTERNAL  BIGA, CKINMI, ERRMSG, MIPRNT, SMALL1, SMALL2, TESTMI
    ..
    .. Intrinsic Functions ..
    
    INTRINSIC ABS, AIMAG, CMPLX, CONJG, COS, MAX, MIN, REAL, SIN
    ..
    SAVE      PASS1
    
    .. Statement Functions ..
    
    REAL      SQ
    ..
    .. Statement Function definitions ..
    
    SQ( CTMP ) = REAL( CTMP )**2 + AIMAG( CTMP )**2
    ..
    DATA      PASS1 /.TRUE./
    """
        
    def SQ( CTMP ):
        np.real( CTMP )**2 + np.imag( CTMP )**2
        
    PASS1 = True
    
        
    #Save some input variables and replace them with values needed to do the self-test
    
    if PASS1: 
        TESTMI( False, XX, CREFIN, MIMCUT, PERFCT, ANYANG, NUMANG, XMU, QEXT,
                QSCA, GQSC, SFORW, SBACK, S1, S2, TFORW, TBACK )

        #PRINT CHECKPOINT 1 OF 8 IN MIEV0
        print "MIEV0 CHECKPOINT 1/8"
            
    #The while loops serve as the go-to/continues labeled 10 and 100 in the Fortran code
    while True:
        
        while True:
            
            #Check input and calculate certain variables from input 
            
            CKINMI( NUMANG, MAXANG, XX, PERFCT, CREFIN, NMOM, ANYANG, XMU )
            
            
            if PERFCT and XX <= 0.1:
            #Use totally-reflecting small-particle limit   
    
                SMALL1( XX, NUMANG, XMU, QEXT, QSCA, GQSC, SFORW, SBACK, 
                        S1, S2, TFORW, TBACK, LITA, LITB ) 
                
                NTRM = 2
                
                #The break serves as the go-to/continue labeled 100 in the Fortran code
                break
                    
                        
            NOABS = True
    
            if not PERFCT:
                
                CIOR = CREFIN
                
                if np.imag( CIOR ) > 0.0:
                    CIOR = np.conjugate( CIOR )
    
                MRE    =   np.real( CIOR )
                MIM    = - np.imag( CIOR )
                NOABS  = MIM <= MIMCUT
                CIORIV = 1.0 / CIOR
                RIORIV = 1.0 / MRE
    
                if XX*max( 1.0, abs( CIOR ) ) <= 0.1:
                    
                    #Use general-refractive-index small-particle limit ( Ref. 2, p. 1508 )
    
                    SMALL2( XX, CIOR, not NOABS, NUMANG, XMU, QEXT, QSCA,
                            GQSC, SFORW, SBACK, S1, S2, TFORW, TBACK, LITA, LITB )
                    
                    NTRM = 2
                    
                    #The break serves as the go-to/continue labeled 100 in the Fortran code
                    break

                
            #PRINT CHECKPOINT 2 OF 8 IN MIEV0
            print "MIEV0 CHECKPOINT 2/8"
                            
            NANGD2 = ( NUMANG + 1 ) / 2
            YESANG = NUMANG > 0
            
            #Estimate number of terms in Mie series ( Ref. 2, p. 1508 )
            if XX <= 8.0:
                
                NTRM = XX + 4.*XX**ONETHR + 1.
    
            elif XX < 4200:
                
                NTRM = XX + 4.05*XX**ONETHR + 2.
    
            else:
                
                NTRM = XX + 4.*XX**ONETHR + 2.
                
    
            if NTRM + 1 > MAXTRM:
                ErrMsg( 'MIEV0--PARAMETER MaxTrm TOO SMALL', True )
    
            #Calculate logarithmic derivatives of J-Bessel-fcn., A-sub-(1 to NTrm)
    
            if not PERFCT:
                BIGA( CIOR, XX, NTRM, NOABS, YESANG, RBIGA, CBIGA )

            #PRINT CHECKPOINT 3 OF 8 IN MIEV0
            print "MIEV0 CHECKPOINT 3/8"
                  
            #Initialize Ricatti-Bessel functions (psi,chi,zeta)-sub-(0,1) for 
            #upward recurrence ( Ref. 1, Eq. 19 )
            XINV   = 1.0 / XX
            PSINM1 = np.sin( XX )
            CHINM1 = np.cos( XX )
            PSIN   = PSINM1*XINV - CHINM1
            CHIN   = CHINM1*XINV + PSINM1
            ZETNM1 = complex( PSINM1, CHINM1 )
            ZETN   = complex( PSIN, CHIN )
            
            #Initialize previous coefficients for -GQSC- series
            ANM1 = 0.0 + 0.0j
            BNM1 = 0.0 + 0.0j
            
            #Initialize angular function  pi and sums for S+, S- ( Ref. 2, p. 1507 )
            if ANYANG:
                
                for J in range( 0, NUMANG ):
                    PINM1[ J ] = 0.0
                    PIN[ J ]   = 1.0
                    SP[ J ]    = 0.0 + 0.0j
                    SM[ J ]    = 0.0 + 0.0j
        
            else:
                
                for J in range( 0, NANGD2 ):
                    PINM1[ J ] = 0.0
                    PIN[ J ]   = 1.0
                    SP [ J ]   = 0.0 + 0.0j
                    SM[ J ]    = 0.0 + 0.0j
                    SPS[ J ]   = 0.0 + 0.0j
                    SMS[ J ]   = 0.0 + 0.0j
          
        
            #Initialize Mie sums for efficiencies, etc.
            QSCA   = 0.0
            GQSC   = 0.0
            SFORW  = 0. + 0.j
            SBACK  = 0. + 0.j
            TFORW[ 0 ] = 0. + 0.j
            TBACK[ 0 ] = 0. + 0.j
        
        
    # ---------  LOOP TO SUM MIE SERIES  -----------------------------------

            #PRINT CHECKPOINT 4 OF 8 IN MIEV0
            print "MIEV0 CHECKPOINT 4/8"
            
            MM     = 1.0
            SPIKE  = 1.0
    
            for N in range( 0, NTRM ):
            #Compute various numerical coefficients
                FN     = N + 1
                RN     = 1.0 / FN
                NP1DN  = 1.0 + RN
                TWONP1 = 2*(N + 1) + 1
                COEFF  = TWONP1 / ( FN*( N + 1 + 1 ) )
                TCOEF  = TWONP1*( FN*( N + 1 + 1 ) )
                
                #Calculate Mie series coefficients
                if PERFCT:
                    #Totally-reflecting case
            
                    AN  = ( ( FN*XINV )*PSIN - PSINM1 ) / ( ( FN*XINV )*ZETN - ZETNM1 )
                    BN  = PSIN / ZETN
    
                elif NOABS:
                    #No-absorption case
    
                    CDENAN = ( RIORIV*RBIGA[ N ] + ( FN*XINV ) ) * ZETN - ZETNM1
                    AN =   ( ( RIORIV*RBIGA[ N ] + ( FN*XINV ) ) * PSIN - PSINM1 ) / CDENAN
                    CDENBN = (  MRE * RBIGA[ N ] + ( FN*XINV ) ) * ZETN - ZETNM1
                    BN =   ( (  MRE * RBIGA[ N ] + ( FN*XINV ) ) * PSIN - PSINM1 ) / CDENBN
    
                else:
                    #Absorptive case
    
                    CDENAN = ( CIORIV * CBIGA[ N ] + (FN*XINV) ) * ZETN - ZETNM1
                    CDENBN = (   CIOR * CBIGA[ N ] + (FN*XINV) ) * ZETN - ZETNM1
                    AN =   ( ( CIORIV * CBIGA[ N ] + (FN*XINV) ) * PSIN - PSINM1 ) / CDENAN
                    BN =   ( (   CIOR * CBIGA[ N ] + (FN*XINV) ) * PSIN - PSINM1 ) / CDENBN
    
                    QSCA = QSCA + TWONP1 * ( SQ( AN ) + SQ( BN ) )
    
                
                if not PERFCT and N > XX:
                    #Flag resonance spikes
                    DENAN  = abs( CDENAN )
                    DENBN  = abs( CDENBN )
                    RATIO  = DENAN / DENBN
    
                    if RATIO <= 0.2 or RATIO >= 5.0:
                        SPIKE = min( SPIKE, DENAN, DENBN )
        
                #Increment Mie sums for non-angle-dependent quantities
    
                SFORW      = SFORW      + TWONP1 * ( AN + BN )
                TFORW[ 0 ] = TFORW[ 0 ] + TCOEF  * ( AN - BN )
                SBACK      = SBACK      + ( MM * TWONP1 ) * ( AN - BN )
                TBACK[ 0 ] = TBACK[ 0 ] + ( MM * TCOEF )  * ( AN + BN )
                GQSC       = GQSC + (FN - RN) * np.real( ANM1 * np.conjugate( AN ) + BNM1 * \
                             np.conjugate( BN ) ) + COEFF * np.real( AN * np.conjugate( BN ) )
                             
                #PRINT CHECKPOINT 5 OF 8 IN MIEV0
                print "MIEV0 CHECKPOINT 5/8"
    
                if YESANG:
                    #Put Mie coefficients in form needed for computing 
                    #S+, S- ( Ref. 2, p. 1507 )
                    ANP = COEFF * ( AN + BN )
                    BNP = COEFF * ( AN - BN )
                    
                    #Increment Mie sums for S+, S- while upward recursing 
                    #angular functions pi and tau
                    if ANYANG:
                        #Arbitrary angles
    
                        #vectorizable loop
                        for J in range( 0, NUMANG ):
                            RTMP = ( XMU[ J ] * PIN[ J ] ) - PINM1[ J ]
                            TAUN =  FN * RTMP - PINM1[ J ]
                            SP[ J ]  = SP[ J ] + ANP * ( PIN[ J ] + TAUN )
                            SM[ J ]  = SM[ J ] + BNP * ( PIN[ J ] - TAUN )
                            PINM1[ J ] = PIN[ J ]
                            PIN[ J ] = ( XMU[ J ] * PIN[ J ] ) + NP1DN * RTMP
    
                    else:
                        #Angles symmetric about 90 degrees
                        ANPM   = MM * ANP
                        BNPM   = MM * BNP
                        #vectorizable loop
                        for J in range( 0, NANGD2 ):
                            RTMP = ( XMU[ J ] * PIN[ J ] ) - PINM1[ J ]
                            TAUN =  FN * RTMP - PINM1[ J ]
                            SP [ J ] = SP [ J ] +  ANP * ( PIN[ J ] + TAUN )
                            SMS[ J ] = SMS[ J ] + BNPM * ( PIN[ J ] + TAUN )
                            SM [ J ] = SM [ J ] +  BNP * ( PIN[ J ] - TAUN )
                            SPS[ J ] = SPS[ J ] + ANPM * ( PIN[ J ] - TAUN )
                            PINM1[ J ] = PIN[ J ]
                            PIN[ J ] = ( XMU[ J ] * PIN[ J ] ) + NP1DN * RTMP
                          
                        
                #Update relevant quantities for next pass through loop
                
                MM     = -MM
                ANM1   = AN
                BNM1   = BN
                
                #Upward recurrence for Ricatti-Bessel functions ( Ref. 1, Eq. 17 )
                
                ZET    = ( TWONP1*XINV )*ZETN - ZETNM1
                ZETNM1 = ZETN
                ZETN   = ZET
                PSINM1 = PSIN
                PSIN   = np.real( ZETN )
    
    #---------- END LOOP TO SUM MIE SERIES --------------------------------

            #PRINT CHECKPOINT 6 OF 8 IN MIEV0
            print "MIEV0 CHECKPOINT 6/8"
    
            QEXT = 2./ XX**2*np.real( SFORW )
    
            if PERFCT or NOABS:
                
                QSCA = QEXT
    
            else:
                
                QSCA = 2./ XX**2 * QSCA
    
                
            GQSC  = 4./ XX**2 * GQSC
            SFORW = 0.5 * SFORW
            SBACK = 0.5 * SBACK
            TFORW[ 1 ] = 0.5 * (   SFORW + 0.25 * TFORW[ 0 ] )
            TFORW[ 0 ] = 0.5 * (   SFORW - 0.25 * TFORW[ 0 ] )
            TBACK[ 1 ] = 0.5 * (   SBACK + 0.25 * TBACK[ 0 ] )
            TBACK[ 0 ] = 0.5 * ( - SBACK + 0.25 * TBACK[ 0 ] )
    
            
            if YESANG:
            #Recover scattering amplitudes from S+, S- ( Ref. 1, Eq. 11 )
            
                if ANYANG:
                    #vectorizable loop
                    for J in range( 0, NUMANG ):
                        S1[ J ] = 0.5*( SP[ J ] + SM[ J ] )
                        S2[ J ] = 0.5*( SP[ J ] - SM[ J ] )
    
                else:
                    #vectorizable loop
                    for J in range( 0, NANGD2 ):
                        S1[ J ] = 0.5*( SP[ J ] + SM[ J ] )
                        S2[ J ] = 0.5*( SP[ J ] - SM[ J ] )
                        
                    #vectorizable loop
                    for J in range( 0, NANGD2 ):
                        S1[ NUMANG - J ] = 0.5*( SPS[ J ] + SMS[ J ] )
                        S2[ NUMANG - J ] = 0.5*( SPS[ J ] - SMS[ J ] )

            break
            #end of while loop representing the go-to/continue labeled 100 in the Fortran code

            
        #PRINT CHECKPOINT 7 OF 8 IN MIEV0
        print "MIEV0 CHECKPOINT 7/8"
                       
        if np.imag( CREFIN ) > 0.0:
            #Take complex conjugates of scattering amplitudes
                
            SFORW  = np.conjugate( SFORW )
            SBACK  = np.conjugate( SBACK )
    
            for I in range( 0, 2 ):
                TFORW[ I ] = np.conjugate( TFORW[ I ] )
                TBACK[ I ] = np.conjugate( TBACK[ I ] )
    
            for J in range( 0, NUMANG ):
                S1[ J ] = np.conjugate( S1[ J ] )
                S2[ J ] = np.conjugate( S2[ J ] )
                    
                            
        if PASS1:
            #Compare test case results with correct answers and abort if bad; 
            #otherwise restore user input and proceed
    
            TESTMI( True, XX, CREFIN, MIMCUT, PERFCT, ANYANG, NUMANG, XMU,
                    QEXT, QSCA, GQSC, SFORW, SBACK, S1, S2, TFORW, TBACK )
                    
            PASS1 = False
            
            
        else:
            if PRNT[ 0 ] or PRNT[ 1 ]:
                MIPRNT( PRNT, XX, PERFCT, CREFIN, NUMANG, XMU, QEXT, QSCA, GQSC, SFORW,
                        SBACK, TFORW, TBACK, S1, S2 )
            break
            #end of while loop representing go-to/continue labeled 10 in Fortran code

    #PRINT CHECKPOINT 8 OF 8 IN MIEV0
    print "MIEV0 CHECKPOINT 8/8"
    
    return
            

  
def CKINMI( NUMANG, MAXANG, XX, PERFCT, CREFIN, NMOM, ANYANG, XMU ):
    
    """     
        Check for bad input to MIEV0
        (NoPMOM version)
    
    Routines called :  ERRMSG, WRTBAD, WRTDIM
        
        
    IMPLICIT NONE
        
    .. Scalar Arguments ..
    
    LOGICAL   ANYANG, PERFCT
    INTEGER   MAXANG, NMOM, NUMANG
    REAL      XX
    COMPLEX   CREFIN
    ..
    .. Array Arguments ..
    
    REAL      XMU( * )
    ..
    .. Local Scalars ..
    
    LOGICAL   INPERR
    INTEGER   I
    ..
    .. External Functions ..
    
    LOGICAL   WRTBAD, WRTDIM
    EXTERNAL  WRTBAD, WRTDIM
    ..
    .. External Subroutines ..
    
    EXTERNAL  ERRMSG
    ..
    .. Intrinsic Functions ..
    
    INTRINSIC ABS, AIMAG, REAL
    ..
    """    
    
        
    INPERR = False
    
    if NUMANG > MAXANG:
        INPERR = WRTDIM( 'MaxAng', NUMANG ) 
    if NUMANG < 0:
        INPERR = WRTBAD( 'NUMANG' )
    
    if XX < 0.:
        INPERR = WRTBAD( 'XX' )
    
    if not PERFCT and np.real( CREFIN ) <= 0.:
        INPERR = WRTBAD( 'CREFIN' )
        
        
    if NMOM != 0:
    
            INPERR = WRTBAD( 'NMOM' )
            print 'For nonzero NMOM, use complete version of MIEVO'
    
        
    if ANYANG:
        #Allow for slight imperfections in computation of cosine
        for I in range( 0, NUMANG ):
    
            if XMU[ I ] < -1.00001 or XMU[ I ] > 1.00001:
                INPERR = WRTBAD( 'XMU' )
    
    else:
            
        for I in range( 0, ( NUMANG + 1 ) / 2 ):
    
            if XMU[ I ] < -0.00001 or XMU[ I ] > 1.00001:
                INPERR = WRTBAD( 'XMU' )
                    
                    
    if INPERR: 
        ErrMsg( 'MIEV0--INPUT ERROR(S).  Aborting...', True )
    
    if XX > 20000.0 or np.real( CREFIN ) > 10.0 or abs( np.imag( CREFIN ) ) > 10.0:
        ErrMsg( 'MIEV0--XX or CREFIN outside tested range', False )

    return



def BIGA( CIOR, XX, NTRM, NOABS, YESANG, RBIGA, CBIGA ):
    
    """
        Calculate logarithmic derivatives of J-Bessel-function
        
    Input :  CIOR, XX, NTRM, NOABS, YESANG  (defined in MIEV0)
    Output :  RBIGA or CBIGA  (defined in MIEV0)
    Routines called :  CONFRA
    
    INTERNAL VARIABLES :
    
        CONFRA     Value of Lentz continued fraction for cBigA(NTrm),
                   used to initialize downward recurrence.
        DOWN       = True, use down-recurrence.  False, do not.
        F1,F2,F3   Arithmetic statement functions used in determining
                   whether to use up-  or down-recurrence ( Ref. 2, Eqs. 6-8 )
        MRE        Real refractive index
        MIM        Imaginary refractive index
        REZINV     1 / ( MRE * XX ); temporary variable for recurrence
        ZINV       1 / ( CIOR * XX ); temporary variable for recurrence
    
    
    IMPLICIT NONE
    
    .. Scalar Arguments ..
    
    LOGICAL   NOABS, YESANG
    INTEGER   NTRM
    REAL      XX
    COMPLEX   CIOR
    ..
    .. Array Arguments ..
    
    REAL      RBIGA( * )
    COMPLEX   CBIGA( * )
    ..
    .. Local Scalars ..
    
    LOGICAL   DOWN
    INTEGER   N
    REAL      MIM, MRE, REZINV, RTMP
    COMPLEX   CTMP, ZINV
    ..
    .. External Functions ..
    
    COMPLEX   CONFRA
    EXTERNAL  CONFRA
    ..
    .. Intrinsic Functions ..
    
    INTRINSIC ABS, AIMAG, COS, EXP, REAL, SIN
    ..
    .. Statement Functions ..
    
    REAL      F1, F2, F3
    ..
    .. Statement Function definitions ..
    
    F1( MRE ) = -8.0 + MRE**2*(26.22 + MRE*(-0.4474 + MRE**3*( 0.00204 - 0.000175*MRE )))
      
    F2( MRE ) = 3.9 + MRE*( -10.8 + 13.78*MRE )
          
    F3( MRE ) = -15.04 + MRE*( 8.42 + 16.35*MRE )
    ..
    """    
    
    def F1( MRE ):
        -8.0 + MRE**2*( 26.22 + MRE*( -0.4474 + MRE**3*( 0.00204 - 0.000175*MRE ) ) )
    def F2( MRE ):
        3.9 + MRE*( -10.8 + 13.78*MRE )
    def F3( MRE ):
        -15.04 + MRE*( 8.42 + 16.35*MRE )
    
    #Decide whether BigA can be calculated by up-recurrence
    MRE = np.real( CIOR )
    MIM = abs( np.imag( CIOR ) )
    
    if MRE < 1.0 or MRE > 10.0 or MIM > 10.0:
            
        DOWN = True
    
    elif YESANG:
            
        DOWN = True
        if MIM*XX < F2( MRE ): 
            DOWN = False
    
    else:
            
        DOWN = True
        if MIM*XX < F1( MRE ):
            DOWN = False
    
                
    ZINV   = 1.0 / ( CIOR*XX )
    REZINV = 1.0 / ( MRE*XX )
    
        
    if DOWN:
    #Compute initial high-order BigA using Lentz method ( Ref. 1, pp. 17-20 )
    
        CTMP = CONFRA( NTRM, ZINV )
    
        #Downward recurrence for BigA ( Ref. 1, Eq. 22 )
        if NOABS:
            #No-absorption case
            RBIGA[ NTRM - 1 ] = np.real( CTMP )
                
            for N in range( NTRM - 1, 2 - 1 - 1, -1 ):
                RBIGA[ N - 1 ] = ( (N + 1)*REZINV ) - 1.0 / ( ( (N + 1)*REZINV ) + RBIGA[ N ] )
    
        else:
            #Absorptive case
            CBIGA[ NTRM - 1 ] = CTMP
    
            for N in range(NTRM - 1, 2 - 1 - 1, -1):
                CBIGA[ N - 1 ] = ( (N + 1)*ZINV ) - 1.0 / ( ( (N + 1)*ZINV ) + CBIGA[ N ] )
    
                    
    else:
        #Upward recurrence for BigA ( Ref. 1, Eqs. 20-21 )
        if NOABS:
            #No-absorption case
            RTMP   = np.sin( MRE*XX )
            RBIGA[ 0 ] = -REZINV + RTMP / ( RTMP*REZINV - np.cos( MRE*XX ) )
    
            for N in range( 2 - 1, NTRM ):
                RBIGA[ N ] = -( (N + 1)*REZINV ) + 1.0 / ( ( (N + 1)*REZINV ) - RBIGA[ N - 1 ] )
    
        else:
            #Absorptive case
            CTMP = np.exp( - complex( 0., 2. ) * CIOR * XX )
            CBIGA[ 0 ] = - ZINV + ( 1.-CTMP ) / ( ZINV * ( 1.-CTMP ) - \
                         complex( 0., 1. )*( 1.+CTMP ) )
    
            for N in range( 2 - 1, NTRM ):
                CBIGA[ N ] = - ( (N + 1)*ZINV ) + 1.0 / ( ( (N + 1)*ZINV ) - CBIGA[ N-1 ] )

    return

  
            
def CONFRA( N, ZINV ):
        
    """       
        Compute Bessel function ratio A-sub-N from its
        continued fraction using Lentz method ( Ref. 1, pp. 17-20 )
    
        ZINV = Reciprocal of argument of A
    
    
    I N T E R N A L    V A R I A B L E S
    ------------------------------------
    
    CAK      Term in continued fraction expansion of A (Ref. 1, Eq. 25)
    CAPT     Factor used in Lentz iteration for A (Ref. 1, Eq. 27)
    CDENOM   Denominator in capT  ( Ref. 1, Eq. 28B )
    CNUMER   Numerator   in capT  ( Ref. 1, Eq. 28A )
    CDTD     Product of two successive denominators of capT factors ( Ref. 1, Eq. 34C )
    CNTN     Product of two successive numerators of capT factors ( Ref. 1, Eq. 34B )
    EPS1     Ill-conditioning criterion
    EPS2     Convergence criterion
    KK       Subscript k of cAk  ( Ref. 1, Eq. 25B )
    KOUNT    Iteration counter ( used only to prevent runaway )
    MAXIT    Max. allowed no. of iterations
    MM       + 1  and - 1, alternately
    
    
    IMPLICIT NONE   
      
    .. Scalar Arguments ..
    
    INTEGER   N
    COMPLEX   ZINV
    ..
    .. Local Scalars ..
    
    INTEGER   KK, KOUNT, MAXIT, MM
    REAL      EPS1, EPS2
    COMPLEX   CAK, CAPT, CDENOM, CDTD, CNTN, CNUMER
    ..
    .. External Subroutines ..
    
    EXTERNAL  ERRMSG
    ..
    .. Intrinsic Functions ..
    
    INTRINSIC ABS, AIMAG, REAL
    ..
    DATA      EPS1 / 1.E-2 / , EPS2 / 1.E-8 /
    DATA      MAXIT / 10000 /
    """    
    
    EPS1 = 1.E-2
    EPS2 = 1.E-8
        
    MAXIT = 10000 
    
        
    #Ref. 1, Eqs. 25a, 27
    CONFRA = ( N + 1 ) * ZINV
    MM     = - 1
    KK     = 2*N + 3
    CAK    = ( MM*KK ) * ZINV
    CDENOM = CAK
    CNUMER = CDENOM + 1.0 / CONFRA
    KOUNT  = 1
    
    #The while loop serves as the go-to/continue labeled 10 in the Fortran code  
    while True:
        KOUNT = KOUNT + 1
            
        if KOUNT > MAXIT:
            ErrMsg( 'ConFra--Iteration failed to converge', True )
    
        #Ref. 2, Eq. 25b
        MM  = - MM
        KK  = KK + 2
        CAK = ( MM*KK ) * ZINV
        #Ref. 2, Eq. 32
        if abs( CNUMER / CAK ) <= EPS1 or abs( CDENOM / CAK ) <= EPS1:
                
        #Ill-conditioned case -- stride two terms instead of one
            
            #Ref. 2, Eqs. 34
            CNTN   = CAK * CNUMER + 1.0
            CDTD   = CAK * CDENOM + 1.0
            CONFRA = ( CNTN / CDTD ) * CONFRA
            #Ref. 2, Eq. 25b
            MM  = - MM
            KK  = KK + 2
            CAK = ( MM*KK ) * ZINV
            #Ref. 2, Eqs. 35
            CNUMER = CAK + CNUMER / CNTN
            CDENOM = CAK + CDENOM / CDTD
            KOUNT  = KOUNT + 1 
                
        else:    
        #Well-conditioned case
            
            #Ref. 2, Eqs. 26, 27
            CAPT   = CNUMER / CDENOM
            CONFRA = CAPT * CONFRA 
            #Check for convergence ( Ref. 2, Eq. 31 )
                
            if abs( np.real( CAPT ) - 1.0 ) >= EPS2 or abs( np.imag( CAPT ) ) >= EPS2:
    
                #Ref. 2, Eqs. 30A-B
                CNUMER = CAK + 1.0 / CNUMER
                CDENOM = CAK + 1.0 / CDENOM
            
            else:
                break

    return

                    
                    
def MIPRNT( PRNT, XX, PERFCT, CREFIN, NUMANG, XMU, QEXT, QSCA, GQSC, SFORW,
            SBACK, TFORW, TBACK, S1, S2 ):
    
    """
        Print scattering quantities of a single particle
        (NoPMOM version)
        
    
    IMPLICIT NONE
    
    .. Scalar Arguments ..
    
    LOGICAL   PERFCT
    INTEGER   NUMANG
    REAL      GQSC, QEXT, QSCA, XX
    COMPLEX   CREFIN, SBACK, SFORW
    ..
    .. Array Arguments ..
    
    LOGICAL   PRNT( * )
    REAL      XMU( * )
    COMPLEX   S1( * ), S2( * ), TBACK( * ), TFORW( * )
    ..
    .. Local Scalars ..
        
    INTEGER   I
    REAL      I1, I2
    ..
    .. Intrinsic Functions ..
    
    INTRINSIC AIMAG, CONJG, REAL
    ..
    """    
    
    if PERFCT: 
        print 'Perfectly Conducting Case, size parameter =', XX
    
    if not PERFCT: 
        print 'Refractive Index:  Real ', np.real( CREFIN ), '  Imag ', np.imag( CREFIN ), \
              ',  Size Parameter =', XX
    
    
    if PRNT[ 0 ] and NUMANG > 0:
    
        print '    cos(angle)  ------- S1 ---------  ------- S2 ---------'
        print '  --- S1*conjg(S2) ---   i1=S1**2   i2=S2**2  (i1+i2)/2'
        print '  DEG POLZN'
            
        for I in range( 1 - 1, NUMANG ):
            I1     = np.real( S1[ I ] )**2 + np.imag( S1[ I ] )**2
            I2     = np.real( S2[ I ] )**2 + np.imag( S2[ I ] )**2
            print  (I + 1), XMU[ I ], S1[ I ], S2[ I ], S1[ I ]*np.conjugate( S2[ I ] ), \
                   I1, I2, 0.5*( I1+I2 ), ( I2-I1 )/( I2+I1 )
                    
                    
    if PRNT[ 1 ]:
    
        print '  Angle ', ' S-sub-1', 'T-sub-1', 'T-sub-2',
        print     0.0 + '  ', SFORW + '  ', TFORW[0] + '  ', TFORW[1] + '  ',
        print    180. + '  ', SBACK + '  ', TBACK[0] + '  ', TBACK[1] + '  '
             
        print ' Efficiency Factors,  extinction:', QEXT
        print '   scattering:', QSCA
        print '   absorption:', QEXT-QSCA
        print '   rad. pressure:', QEXT-GQSC

    return

                    
                    
def SMALL1( XX, NUMANG, XMU, QEXT, QSCA, GQSC, SFORW, SBACK, S1, S2, TFORW, TBACK, A, B ):
    
    """
        Small-particle limit of Mie quantities in totally reflecting
        limit ( Mie series truncated after 2 terms )
    
        A,B        First two Mie coefficients, with numerator and
                   denominator expanded in powers of XX ( a factor
                   of XX**3 is missing but is restored before return
                   to calling program )  ( Ref. 2, p. 1508 )
    
    
    IMPLICIT NONE
    
    .. Parameters ..
    
    REAL      TWOTHR, FIVTHR, FIVNIN
    PARAMETER ( TWOTHR = 2. / 3., FIVTHR = 5. / 3., FIVNIN = 5. / 9. )
    ..
    .. Scalar Arguments ..
    
    INTEGER   NUMANG
    REAL      GQSC, QEXT, QSCA, XX
    COMPLEX   SBACK, SFORW
    ..
    .. Array Arguments ..
    
    REAL      XMU( * )
    COMPLEX   A( * ), B( * ), S1( * ), S2( * ), TBACK( * ), TFORW( * )
    ..
    .. Local Scalars ..
    
    INTEGER   J
    REAL      RTMP
    COMPLEX   CTMP
    ..
    .. Intrinsic Functions ..
    
    INTRINSIC AIMAG, CMPLX, CONJG, REAL
    ..
    .. Statement Functions ..
    
    REAL      SQ
    ..
    .. Statement Function definitions ..
    
    SQ( CTMP ) = REAL( CTMP )**2 + AIMAG( CTMP )**2
    ..
    """  
        
    def SQ( CTMP ):
        np.real( CTMP )**2 + np.imag( CTMP )**2
    
    
    A[ 0 ] = complex( 0., TWOTHR*( 1. - 0.2*XX**2 ) ) / complex( 1. - 0.5*XX**2, TWOTHR*XX**3 )
    B[ 0 ] = complex( 0., -( 1. - 0.1*XX**2 ) / 3. ) / complex( 1. + 0.5*XX**2, -XX**3 / 3. )
    
    A[ 1 ] = complex( 0.,   XX**2 / 30. )
    B[ 1 ] = complex( 0., - XX**2 / 45. )
    
    QSCA   = 6.* XX**4 *( SQ( A[ 0 ] ) + SQ( B[ 0 ] ) + FIVTHR*( SQ( A[ 1 ] ) + SQ( B[ 1 ] ) ) )
    QEXT   = QSCA
    GQSC   = 6.* XX**4 * np.real( A[ 0 ]*np.conjugate( A[ 1 ] + B[ 0 ] ) + \
             ( B[ 0 ] + FIVNIN*A[ 1 ] )*np.conjugate( B[ 1 ] ) )
    
    RTMP   = 1.5 * XX**3
    SFORW  = RTMP * ( A[ 0 ] + B[ 0 ] + FIVTHR * ( A[ 1 ] + B[ 1 ] ) )
    SBACK  = RTMP * ( A[ 0 ] - B[ 0 ] - FIVTHR * ( A[ 1 ] - B[ 1 ] ) )
    TFORW[ 0 ] = RTMP*( B[ 0 ] + FIVTHR * ( 2.*B[ 1 ] - A[ 1 ] ) )
    TFORW[ 1 ] = RTMP*( A[ 0 ] + FIVTHR * ( 2.*A[ 1 ] - B[ 1 ] ) )
    TBACK[ 0 ] = RTMP*( B[ 0 ] - FIVTHR * ( 2.*B[ 1 ] + A[ 1 ] ) )
    TBACK[ 1 ] = RTMP*( A[ 0 ] - FIVTHR * ( 2.*A[ 1 ] + B[ 1 ] ) )
        
        
    for J in range( 0, NUMANG ):
        S1[ J ] = RTMP*( A[ 0 ] + B[ 0 ]*XMU[ J ] + FIVTHR*( A[ 1 ]*XMU[ J ] + \
                  B[ 1 ]*( 2.*XMU[ J ]**2 - 1.) ) )
        S2[ J ] = RTMP*( B[ 0 ] + A[ 0 ]*XMU[ J ] + FIVTHR*( B[ 1 ]*XMU[ J ] + \
                  A[ 1 ]*( 2.*XMU[ J ]**2 - 1.) ) )
    
    #Recover actual Mie coefficients
    A[ 0 ] = XX**3 * A[ 0 ]
    A[ 1 ] = XX**3 * A[ 1 ]
    B[ 0 ] = XX**3 * B[ 0 ]
    B[ 1 ] = XX**3 * B[ 1 ]
             
    return

    
        
def SMALL2( XX, CIOR, CALCQE, NUMANG, XMU, QEXT, QSCA, GQSC, SFORW,
            SBACK, S1, S2, TFORW, TBACK, A, B ):
        
    """    
        Small-particle limit of Mie quantities for general refractive
        index ( Mie series truncated after 2 terms )
    
        A,B        First two Mie coefficients, with numerator and
                   denominator expanded in powers of XX ( a factor
                   of XX**3 is missing but is restored before return
                   to calling program )  ( Ref. 2, p. 1508 )
    
        CIORSQ     Square of refractive index
    
    
    IMPLICIT NONE
    
    .. Parameters ..
    
    REAL      TWOTHR, FIVTHR
    PARAMETER  ( TWOTHR = 2./3., FIVTHR = 5./3. )
    ..
    .. Scalar Arguments ..
    
    LOGICAL   CALCQE
    INTEGER   NUMANG
    REAL      GQSC, QEXT, QSCA, XX
    COMPLEX   CIOR, SBACK, SFORW
    ..
    .. Array Arguments ..
    
    REAL      XMU( * )
    COMPLEX   A( * ), B( * ), S1( * ), S2( * ), TBACK( * ), TFORW( * )
    ..
    .. Local Scalars ..
    
    INTEGER   J
    REAL      RTMP
    COMPLEX   CIORSQ, CTMP
    ..
    .. Intrinsic Functions ..
    
    INTRINSIC AIMAG, CMPLX, CONJG, REAL
    ..
    .. Statement Functions ..
    
    REAL      SQ
    ..
    .. Statement Function definitions ..
    
    SQ( CTMP ) = REAL( CTMP )**2 + AIMAG( CTMP )**2
    ..
    """   
    
    def SQ( CTMP ):
        np.real( CTMP )**2 + np.imag( CTMP )**2
    
        
    CIORSQ = CIOR**2
    CTMP   = complex( 0., TWOTHR )*( CIORSQ - 1. )
    A[ 0 ] = CTMP*( 1. - 0.1*XX**2 + ( CIORSQ / 350.+ 1./ 280. )*XX**4 ) / \
             ( CIORSQ + 2.+ ( 1.- 0.7*CIORSQ )*XX**2 - ( CIORSQ**2 / 175. - \
             0.275*CIORSQ + 0.25 )*XX**4 + XX**3 * CTMP*( 1.- 0.1*XX**2 ) )
    
    B[ 0 ] = ( XX**2 / 30. ) * CTMP * ( 1.+ ( CIORSQ / 35.- 1./ 14. )*XX**2 ) / \
             ( 1.- ( CIORSQ / 15.- 1./ 6. )*XX**2 )
    
    A[ 1 ] = ( 0.1*XX**2 )*CTMP*( 1.- XX**2 / 14. ) / ( 2.*CIORSQ + 3. - \
             ( CIORSQ / 7.- 0.5 )*XX**2 )
    
    QSCA   = 6.* XX**4 * ( SQ( A[ 0 ] ) + SQ( B[ 0 ] ) + FIVTHR*SQ( A[ 1 ] ) )
    GQSC   = 6.* XX**4 * np.real( A[ 0 ]*np.conjugate( A[ 1 ] + B[ 0 ] ) )
    QEXT   = QSCA
        
    if CALCQE: 
        QEXT  = 6.*XX*np.real( A[ 0 ] + B[ 0 ] + FIVTHR*A[ 1 ] )
    
    RTMP   = 1.5 * XX**3
    SFORW  = RTMP*( A[ 0 ] + B[ 0 ] + FIVTHR*A[ 1 ] )
    SBACK  = RTMP*( A[ 0 ] - B[ 0 ] - FIVTHR*A[ 1 ] )
    TFORW[ 0 ] = RTMP*( B[ 0 ] - FIVTHR*A[ 1 ] )
    TFORW[ 1 ] = RTMP*( A[ 0 ] + 2.*FIVTHR*A[ 1 ] )
    TBACK[ 0 ] = TFORW[ 0 ]
    TBACK[ 1 ] = RTMP*( A[ 0 ] - 2.*FIVTHR*A[ 1 ] )
        
        
    for J in range( 0, NUMANG ):
        S1[ J ] = RTMP*( A[ 0 ] + ( B[ 0 ] + FIVTHR*A[ 1 ] )*XMU[ J ] )
        S2[ J ] = RTMP*( B[ 0 ] + A[ 0 ]*XMU[ J ] + FIVTHR * A[ 1 ]*( 2.*XMU[ J ]**2 - 1. ) )
        
    #Recover actual Mie coefficients
    A[ 0 ] = XX**3 * A[ 0 ]
    A[ 1 ] = XX**3 * A[ 1 ]
    B[ 0 ] = XX**3 * B[ 0 ]
    B[ 1 ] = 0. + 0.j
    
    return


    
def TESTMI( COMPAR, XX, CREFIN, MIMCUT, PERFCT, ANYANG, NUMANG, XMU, QEXT, 
            QSCA, GQSC, SFORW, SBACK, S1, S2, TFORW, TBACK ):
    
    """
        Set up to run test case when  COMPAR = False;  when  = True,
        compare Mie code test case results with correct answers
        and abort if even one result is inaccurate.
        (NoPMOM version)
    
        The test case is :  Mie size parameter = 10
                            refractive index   = 1.5 - 0.1 i
                            scattering angle = 140 degrees
                            1 Sekera moment
    
        Results for this case may be found among the test cases
        at the end of reference (1).
    
        *** NOTE *** When running on some computers, esp. in single
        precision, the Accur criterion below may have to be relaxed.
        However, if Accur must be set larger than 10**-3 for some
        size parameters, your computer is probably not accurate
        enough to do Mie computations for those size parameters.
    
    Routines called :  ERRMSG, MIPRNT, TSTBAD
        
        
    IMPLICIT NONE    
        
    .. Scalar Arguments ..
    
    LOGICAL   ANYANG, COMPAR, PERFCT
    INTEGER   NUMANG
    REAL      GQSC, MIMCUT, QEXT, QSCA, XX
    COMPLEX   CREFIN, SBACK, SFORW
    ..
    .. Array Arguments ..
    
    REAL      XMU( * )
    COMPLEX   S1( * ), S2( * ), TBACK( * ), TFORW( * )
    ..
    .. Local Scalars ..
    
    LOGICAL   ANYSAV, OK, PERSAV
    INTEGER   N, NUMSAV
    REAL      ACCUR, CALC, EXACT, MIMSAV, TESTGQ, TESTQE, TESTQS, XMUSAV, XXSAV
    COMPLEX   CRESAV, TESTS1, TESTS2, TESTSB, TESTSF
    ..
    .. Local Arrays ..
    
    LOGICAL   PRNT( 2 )
    COMPLEX   TESTTB( 2 ), TESTTF( 2 )
    ..
    .. External Functions ..
    
    LOGICAL   TSTBAD
    EXTERNAL  TSTBAD
    ..
    .. External Subroutines ..
    
    EXTERNAL  ERRMSG, MIPRNT
    ..
    .. Intrinsic Functions ..
    
    INTRINSIC ABS, AIMAG, REAL
    ..
    .. Statement Functions ..
    
    LOGICAL WRONG
    ..
    SAVE    XXSAV, CRESAV, MIMSAV, PERSAV, ANYSAV, NUMSAV, XMUSAV
          
    DATA    TESTQE / 2.459791 /,
            TESTQS / 1.235144 /,
            TESTGQ / 1.139235 /,
            TESTSF / ( 61.49476, -3.177994 ) /,
            TESTSB / ( 1.493434, 0.2963657 ) /,
            TESTS1 / ( -0.1548380, -1.128972) /,
            TESTS2 / ( 0.05669755, 0.5425681) /,
            TESTTF / ( 12.95238, -136.6436 ), ( 48.54238, 133.4656 ) /,
            TESTTB / ( 41.88414, -15.57833 ), ( 43.37758, -15.28196 )/
         
    DATA    ACCUR / 1.E-4 /
    ..
    .. Statement Function definitions ..
    
    WRONG( CALC, EXACT ) = ABS( ( CALC - EXACT ) / EXACT ).GT.ACCUR
    ..
    """
    
    TESTQE = 2.459791 
    TESTQS = 1.235144 
    TESTGQ = 1.139235 
    TESTSF = 61.49476 - 3.177994j 
    TESTSB = 1.493434 + 0.2963657j 
    TESTS1 = -0.1548380 - 1.128972j 
    TESTS2 = 0.05669755 + 0.5425681j 
    TESTTF = [ 12.95238 - 136.6436j, 48.54238 + 133.4656j ] 
    TESTTB = [ 41.88414 - 15.57833j, 43.37758 - 15.28196j ]
         
    ACCUR = 1.E-4 
    
    def WRONG( CALC, EXACT ):
        abs( ( CALC - EXACT ) / EXACT ) > ACCUR
    
    
    if not COMPAR:
        #Save certain user input values
        XXSAV  = XX
        CRESAV = CREFIN
        MIMSAV = MIMCUT
        PERSAV = PERFCT
        ANYSAV = ANYANG
        NUMSAV = NUMANG
        XMUSAV = XMU[ 0 ]
        #Reset input values for test case
        XX     = 10.0
        CREFIN = 1.5 - 0.1j
        MIMCUT = 0.0
        PERFCT = False
        ANYANG = True
        NUMANG = 1
        XMU[ 0 ] = -0.7660444
    
    else:
        #Compare test case results with correct answers and abort if bad
        OK     =  True
    
        if WRONG( QEXT, TESTQE ):
            OK =  TSTBAD( 'QEXT', abs( ( QEXT - TESTQE ) / TESTQE ) )
         
        if WRONG( QSCA, TESTQS ):
            OK =  TSTBAD( 'QSCA', abs( ( QSCA - TESTQS ) / TESTQS ) )
         
        if WRONG( GQSC, TESTGQ ):
            OK =  TSTBAD( 'GQSC', abs( ( GQSC - TESTGQ ) / TESTGQ ) )
    
        if WRONG( np.real( SFORW ), np.real( TESTSF ) ) or WRONG( np.imag( SFORW ),
                  np.imag( TESTSF ) ):
            OK = TSTBAD( 'SFORW', abs( ( SFORW - TESTSF ) / TESTSF ) )
    
        if WRONG( np.real( SBACK ), np.real( TESTSB ) ) or WRONG( np.imag( SBACK ),
                  np.imag( TESTSB ) ):
            OK = TSTBAD( 'SBACK', abs( ( SBACK - TESTSB ) / TESTSB ) )
    
        if WRONG( np.real( S1[ 0 ] ), np.real( TESTS1 ) ) or WRONG( np.imag( S1[ 0 ] ),
                  np.imag( TESTS1 ) ):
            OK = TSTBAD( 'S1', abs( ( S1[ 0 ] - TESTS1 ) / TESTS1 ) )
    
        if WRONG( np.real( S2[ 0 ] ), np.real( TESTS2 ) ) or WRONG( np.imag( S2[ 0 ] ),
                  np.imag( TESTS2 ) ):
            OK = TSTBAD( 'S2', abs( ( S2[ 0 ] - TESTS2 ) / TESTS2 ) )
                
                
        for N in range( 0, 2 ):
    
            if WRONG( np.real( TFORW[ N ] ), np.real( TESTTF[ N ] ) ) or \
               WRONG( np.imag( TFORW[ N ] ), np.imag( TESTTF[ N ] ) ):
                OK =  TSTBAD( 'TFORW', abs( ( TFORW[ N ] - TESTTF[ N ] ) / TESTTF[ N ] ) )
                    
            if WRONG( np.real( TBACK[ N ] ), np.real( TESTTB[ N ] ) ) or \
               WRONG( np.imag( TBACK[ N ] ), np.imag( TESTTB[ N ] ) ):
                OK =  TSTBAD( 'TBACK', abs( ( TBACK[ N ] - TESTTB[ N ] ) / TESTTB[ N ] ) )

                    
        if not OK:
                
            PRNT[ 0 ] = True
            PRNT[ 1 ] = True
    
            MIPRNT( PRNT, XX, PERFCT, CREFIN, NUMANG, XMU, QEXT, QSCA,
                           GQSC, SFORW, SBACK, TFORW, TBACK, S1, S2 )
    
            ErrMsg( 'MIEV0 -- Self-test failed', True )
    
        #Restore user input values
        XX     = XXSAV
        CREFIN = CRESAV
        MIMCUT = MIMSAV
        PERFCT = PERSAV
        ANYANG = ANYSAV
        NUMANG = NUMSAV
        XMU[ 0 ] = XMUSAV
        
    return
