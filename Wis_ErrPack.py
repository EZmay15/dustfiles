

#The original Wiscombe ErrPack Fortran file is also in the dustfiles repository for
#reference/comparison as ErrPackcopy.f! I found the code on the page "Mie Scattering" by
#Scott Prahl, who had zipped a version of Wiscombe's code and made it
#downloadable (http://omlc.org/software/mie/)

def ErrMsg( MESSAG, FATAL ):
    
    """
    Print out a warning or error message; abort if error 
    after making symbolic dump (machine-specific)

    .. Scalar Arguments ..

    CHARACTER MESSAG*(*)
    LOGICAL   FATAL
    ..
    .. Local Scalars ..

    LOGICAL   MSGLIM
    INTEGER   MAXMSG, NUMMSG
    ..
    .. External Subroutines ..

    EXTERNAL  SYMDUMP
    ..
    SAVE      MAXMSG, NUMMSG, MSGLIM
    DATA      NUMMSG / 0 /,  MAXMSG / 100 /,  MSGLIM /.FALSE./
    """
    
    NUMMSG = 0  
    MAXMSG = 100 
    MSGLIM = False

    if FATAL:
        print '******ERROR***** ', MESSAG
        
    #Example symbolic dump call for Cray
    #CALL SYMDUMP( '-B -c3' )
    
    NUMMSG = NUMMSG + 1
    
    #while True loop allows MSGLIM condition to end the function with a break
    while True:
        if MSGLIM: 
            break
    
        if NUMMSG <= MAXMSG:
            print '******WARNING***** ', MESSAG
            break
        
        else:
            print '****** TOO MANY WARNING MESSAGES -- They will no longer be printed *******'
            MSGLIM = True
            break



def WrtBad( VarNam ):
    
    """
    Write names of erroneous variables and return 'TRUE'
    INPUT : VarNam = Name of erroneous variable to be written (CHARACTER, any length)

    .. Scalar Arguments ..

    CHARACTER VarNam*(*)
    ..
    .. Local Scalars ..

    INTEGER   MAXMSG, NUMMSG
    ..
    .. External Subroutines ..

    EXTERNAL  ErrMsg
    ..
    SAVE      NUMMSG, MAXMSG
    DATA      NUMMSG / 0 /, MAXMSG / 50 /
    """

    NUMMSG = 0 
    MAXMSG = 50 

    NUMMSG = NUMMSG + 1
    print '**** Input variable ', VarNam, ' in error  ****'
    return True

    if NUMMSG == MAXMSG :
        ErrMsg('Too many input errors. Aborting...', True)



def WrtDim( DimNam, Minval ): 
    
    """
    Write name of too-small symbolic dimension and 
    the value it should be increased to; return 'TRUE'

    INPUT: DimNam = string; Name of symbolic dimension which is too small (CHARACTER, any length)
            Minval = int; Value to which that dimension should be increased (at least)

    .. Scalar Arguments ..

    CHARACTER DimNam*(*)
    INTEGER   Minval
    ..
    """
    
    print '**** Symbolic dimension ', DimNam, ' should be increased to at least ', Minval
    return True



def TstBad( VarNam, RelErr ):

    """
    Write name (VarNam) of variable failing self-test and its
    percent error from the correct value;  return  'FALSE'.

    .. Scalar Arguments ..

    CHARACTER VarNam*(*)
    REAL      RelErr
    ..
    """
    
    print '*** Output variable ', VarNam, ' differed by ', RelErr, 
    ' percent from correct value.  Self-test failed.'
    return False

