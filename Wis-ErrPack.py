
# coding: utf-8

# In[10]:

def ErrMsg( MESSAG, FATAL ):
    
    """
c Print out a warning or error message; abort if error 
c after making symbolic dump (machine-specific)

c     .. Scalar Arguments ..

      CHARACTER MESSAG*(*)
      LOGICAL   FATAL
c     ..
c     .. Local Scalars ..

      LOGICAL   MSGLIM
      INTEGER   MAXMSG, NUMMSG
c     ..
c     .. External Subroutines ..

ccccc EXTERNAL  SYMDUMP
c     ..
      SAVE      MAXMSG, NUMMSG, MSGLIM
      DATA      NUMMSG / 0 /,  MAXMSG / 100 /,  MSGLIM /.FALSE./
    """
    
    NUMMSG = 0  
    MAXMSG = 100 
    MSGLIM = False

    if FATAL :
        print '******ERROR***** ', MESSAG
        
#Example symbolic dump call for Cray
#CALL SYMDUMP( '-B -c3' )
    
    NUMMSG = NUMMSG + 1
    
    if MSGLIM : 
        return
    
    if NUMMSG <= MAXMSG :
        print '******WARNING***** ', MESSAG    
    else :
        print '****** TOO MANY WARNING MESSAGES -- They will no longer be printed *******'
        MSGLIM = True


# In[11]:

def WrtBad( VarNam ):
    
    """
c Write names of erroneous variables and return 'TRUE'
c INPUT : VarNam = Name of erroneous variable to be written (CHARACTER, any length)

c     .. Scalar Arguments ..

      CHARACTER VarNam*(*)
c     ..
c     .. Local Scalars ..

      INTEGER   MAXMSG, NUMMSG
c     ..
c     .. External Subroutines ..

      EXTERNAL  ErrMsg
c     ..
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


# In[12]:

def WrtDim( DimNam, Minval ): 
    
    """
c Write name of too-small symbolic dimension and 
c the value it should be increased to; return 'TRUE'

c INPUT: DimNam = string; Name of symbolic dimension which is too small (CHARACTER, any length)
c        Minval = int; Value to which that dimension should be increased (at least)

c     .. Scalar Arguments ..

      CHARACTER DimNam*(*)
      INTEGER   Minval
c     ..
    """
    
    print '**** Symbolic dimension ', DimNam, ' should be increased to at least ', Minval
    return True


# In[13]:

def TstBad( VarNam, RelErr ):

    """
c Write name (VarNam) of variable failing self-test and its
c percent error from the correct value;  return  'FALSE'.

c     .. Scalar Arguments ..

      CHARACTER VarNam*(*)
      REAL      RelErr
c     ..
    """
    
    print '*** Output variable ', VarNam, ' differed by ', RelErr, 
    ' percent from correct value.  Self-test failed.'
    return False

