      SUBROUTINE ENVIRN(OUTCH,MCHEPS,DWARF)
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C  The RK suite requires some environmental parameters that are provided by
C  this subroutine.  The values provided with the distribution codes are those
C  appropriate to the IEEE standard.  They must be altered, if necessary, to
C  those appropriate to the computing system you are using before calling the 
C  codes of the suite.  
C
C        ================================================================
C        ================================================================
C        TO MAKE SURE THAT THESE MACHINE AND INSTALLATION DEPENDENT 
C        QUANTITIES ARE SPECIFIED PROPERLY, THE DISTRIBUTION VERSION 
C        WRITES A MESSAGE ABOUT THE MATTER TO THE STANDARD OUTPUT CHANNEL
C        AND TERMINATES THE RUN.  THE VALUES PROVIDED IN THE DISTRIBUTION
C        VERSION SHOULD BE ALTERED, IF NECESSARY, AND THE "WRITE" AND 
C        "STOP" STATEMENTS COMMENTED OUT.
C        ================================================================
C        ================================================================
C
C  OUTPUT VARIABLES 
C
C     OUTCH     - INTEGER
C                 Standard output channel
C     MCHEPS    - DOUBLE PRECISION
C                 MCHEPS is the largest positive number such that
C                 1.0D0 + MCHEPS = 1.0D0. 
C     DWARF     - DOUBLE PRECISION
C                 DWARF is the smallest positive number.
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C     .. Scalar Arguments ..  
      INTEGER           OUTCH
      DOUBLE PRECISION  DWARF, MCHEPS
C     .. Executable Statements ..      
C
C  The following six statements are to be Commented out after verification that
C  the machine and installation dependent quantities are specified correctly.
C  If you pass copies of RKSUITE on to others, please give them the whole
C  distribution version of RKSUITE, and in particular, give them a version 
C  of ENVIRN that does not have the following six statements Commented out.
c      WRITE(*,*) ' Before using RKSUITE, you must verify that the  '
c      WRITE(*,*) ' machine- and installation-dependent quantities  '
c      WRITE(*,*) ' specified in the subroutine ENVIRN are correct, '
c      WRITE(*,*) ' and then Comment these WRITE statements and the '
c      WRITE(*,*) ' STOP statement out of ENVIRN.                   '
c      STOP
C
C  The following values are appropriate to IEEE arithmetic with the typical
C  standard output channel.
C
      OUTCH = 6
      MCHEPS = 1.11D-16
      DWARF = 2.23D-308
C      
C------------------------------------------------------------------------------
C  If you have the routines D1MACH and I1MACH on your system, you could
C  replace the preceding statements by the following ones to obtain the 
C  appropriate machine dependent numbers. The routines D1MACH and I1MACH 
C  are public domain software.  They are available from NETLIB.
C      .. Scalar Arguments ..  
C      INTEGER           OUTCH
C      DOUBLE PRECISION  DWARF, MCHEPS
C      .. External Functions ..
C      INTEGER           I1MACH
C      DOUBLE PRECISION  D1MACH
C      .. Executable Statements ..
C
C      OUTCH = I1MACH(2)
C      MCHEPS = D1MACH(3)
C      DWARF = D1MACH(1)
C
C  If you have the NAG Fortran Library available on your system, you could 
C  replace the preceding statements by the following ones to obtain the 
C  appropriate machine dependent numbers.
C
C      .. Scalar Arguments ..  
C      INTEGER           OUTCH
C      DOUBLE PRECISION  DWARF, MCHEPS
C      .. External Functions ..
C      DOUBLE PRECISION  X02AJF, X02AMF
C      .. Executable Statements ..
C
C      CALL X04AAF(0,OUTCH)
C      MCHEPS = X02AJF()
C      DWARF = X02AMF()
C
C  If you have the IMSL MATH/LIBRARY available on your system, you could
C  replace the preceding statements by the following ones to obtain the
C  appropriate machine dependent numbers.
C
C      .. Scalar Arguments ..  
C      INTEGER           OUTCH
C      DOUBLE PRECISION  DWARF, MCHEPS
C      .. External Functions ..
C      DOUBLE PRECISION  DMACH
C      .. Executable Statements ..
C
C      CALL UMACH(2,OUTCH)
C      MCHEPS = DMACH(4)
C      DWARF = DMACH(1)
C------------------------------------------------------------------------------
C
      RETURN
      END
