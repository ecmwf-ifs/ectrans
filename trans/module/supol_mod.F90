MODULE SUPOL_MOD
CONTAINS
SUBROUTINE SUPOL(KNSMAX,DDMU,PFN,DDPOL,DDA,DDC,DDD,DDE,DDH,DDI)

!**** *SUPOL * - Routine to compute the Legendre polynomials

!     Purpose.
!     --------
!           For a given value of mu, computes the Legendre
!           polynomials.

!**   Interface.
!     ----------
!        *CALL* *SUPOL(KNSMAX,DDMU,PFN,DDPOL,
!        DDA,DDC,DDD,DDE,DDH,DDI)

!        Explicit arguments :
!        --------------------
!              PFN      :  Fourier coefficients of series expansion 
!                          for the ordinary Legendre polynomials
!              KNSMAX   :  Truncation  (triangular)
!              DDMU     :  Abscissa at which the polynomials are computed (mu)
!              DDPOL    :  Polynomials (the first index is m and the second n)


!        Implicit arguments :   None
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 87-10-15
!        K. YESSAD (MAY 1998): modification to avoid underflow.
!        R. El Khatib 11-Apr-2007 Emulation of vectorized quadruple precision
!                                 on NEC
!        Nils Wedi + Mats Hamrud, 2009-02-05 revised following Swarztrauber, 2002
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KNSMAX
REAL(KIND=JPRB)   ,INTENT(IN)  :: DDMU
REAL(KIND=JPRB)   ,INTENT(IN)  :: PFN(0:KNSMAX,0:KNSMAX)
REAL(KIND=JPRB)   ,INTENT(IN)  :: DDC(0:KNSMAX,0:KNSMAX)
REAL(KIND=JPRB)   ,INTENT(IN)  :: DDD(0:KNSMAX,0:KNSMAX)
REAL(KIND=JPRB)   ,INTENT(IN)  :: DDE(0:KNSMAX,0:KNSMAX)
REAL(KIND=JPRB)   ,INTENT(IN)  :: DDA(0:KNSMAX),DDI(0:KNSMAX),DDH(0:KNSMAX)
REAL(KIND=JPRB)   ,INTENT(OUT) :: DDPOL(0:KNSMAX,0:KNSMAX)

REAL(KIND=JPRB) :: DLX,DLX1,DLSITA,DL1SITA,DLK,DL1,DLS,DLLDN

INTEGER(KIND=JPIM), PARAMETER :: JPKD=KIND(DLX)

INTEGER(KIND=JPIM) :: JM, JN, JK
REAL(KIND=JPRB) :: Z

!     ------------------------------------------------------------------

!*       1. First two columns.
!           ------------------

DLX=DDMU
DLX1=ACOS(DLX)
DLSITA=SQRT(1.0_JPRB-DLX*DLX)

DDPOL(0,0)=1._JPRB
DLLDN = 0.0_JPRB

! IF WE ARE LESS THAN 1Meter FROM THE POLE,
IF(ABS(REAL(DLSITA,KIND(Z))) <= SQRT(EPSILON(Z)))THEN
  DLX=1._JPRB
  DLSITA=0._JPRB
  DL1SITA=0._JPRB
ELSE
  DL1SITA=1.0_JPRB/DLSITA
ENDIF

!*          ordinary Legendre polynomials from series expansion
!           ---------------------------------------------------

! even N
DO JN=2,KNSMAX,2
  DLK = 0.5_JPRB*PFN(JN,0)
  DLLDN = 0.0_JPRB
  ! represented by only even k
  DO JK=2,JN,2
    ! normalised ordinary Legendre polynomial == \overbar{P_n}^0
    DLK = DLK + PFN(JN,JK)*COS(DDI(JK)*DLX1)
    ! normalised associated Legendre polynomial == \overbar{P_n}^1
    DLLDN = DLLDN + DDA(JN)*PFN(JN,JK)*DDI(JK)*SIN(DDI(JK)*DLX1)
  ENDDO
  DDPOL(0,JN) = DLK
  DDPOL(1,JN) = DLLDN
ENDDO
! odd N
DO JN=1,KNSMAX,2
  DLK = 0.0_JPRB
  DLLDN = 0.0_JPRB
  ! represented by only odd k
  DO JK=1,JN,2
    ! normalised ordinary Legendre polynomial == \overbar{P_n}^0
    DLK = DLK + PFN(JN,JK)*COS(DDI(JK)*DLX1)
    ! normalised associated Legendre polynomial == \overbar{P_n}^1
    DLLDN = DLLDN + DDA(JN)*PFN(JN,JK)*DDI(JK)*SIN(DDI(JK)*DLX1)
  ENDDO
  DDPOL(0,JN) = DLK
  DDPOL(1,JN) = DLLDN
ENDDO

!     ------------------------------------------------------------------

!*       2. Diagonal (the terms 0,0 and 1,1 have already been computed)
!           Belousov, equation (23)
!           -----------------------------------------------------------

DLS=DL1SITA*TINY(DLS)

#ifdef VPP
!OCL SCALAR
#endif
DO JN=2,KNSMAX
  DDPOL(JN,JN)=DDPOL(JN-1,JN-1)*DLSITA*DDH(JN)
  IF ( ABS(DDPOL(JN,JN))  <  DLS ) DDPOL(JN,JN)=0.0_JPRB
ENDDO

!     ------------------------------------------------------------------

!*       3. General recurrence (Belousov, equation 17)
!           -----------------------------------------

DO JN=3,KNSMAX
!DIR$ IVDEP
!OCL NOVREC
  DO JM=2,JN-1
    DDPOL(JM,JN)=DDC(JM,JN)*DDPOL(JM-2,JN-2)&
     &-DDD(JM,JN)*DDPOL(JM-2,JN-1)*DLX &
     &+DDE(JM,JN)*DDPOL(JM  ,JN-1)*DLX
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE SUPOL
END MODULE SUPOL_MOD
