MODULE SUPOL_MOD
CONTAINS
SUBROUTINE SUPOL(KNSMAX,PDDMU,PFN,PDDPOL,PDDA,PDDC,PDDD,PDDE,PDDH,PDDI)

!**** *SUPOL * - Routine to compute the Legendre polynomials

!     Purpose.
!     --------
!           For a given value of mu, computes the Legendre polynomials.

!**   Interface.
!     ----------
!        *CALL* *SUPOL(...)

!        Explicit arguments :
!        --------------------
!        KNSMAX    :  Truncation  (triangular)                            [in]
!        PDDMU     :  Abscissa at which the polynomials are computed (mu) [in]
!        PFN       :  Fourier coefficients of series expansion 
!                     for the ordinary Legendre polynomials               [in]
!        PDDPOL    :  Polynomials (the first index is m and the second n) [out]
!        PDDA to PDDI: intermediate precomputed quantities (see caller)   [in]

!        Implicit arguments :   None
!        --------------------

!     Method.
!     -------
!        See documentation about spectral transforms
!         (doc (IDTS) by K. Yessad, appendix 3, or doc (NTA30) by M. Rochas)

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
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        R. El Khatib 11-Apr-2007 Emulation of vectorized quadruple precision
!                                 on NEC
!        K. YESSAD (NOV 2008): make consistent arp/SUPOLA and tfl/SUPOL.
!        Nils Wedi + Mats Hamrud, 2009-02-05 revised following Swarztrauber, 2002
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KNSMAX
REAL(KIND=JPRB)   ,INTENT(IN)  :: PDDMU
REAL(KIND=JPRB)   ,INTENT(IN)  :: PFN(0:KNSMAX,0:KNSMAX)
REAL(KIND=JPRB)   ,INTENT(OUT) :: PDDPOL(0:KNSMAX,0:KNSMAX)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PDDA(0:KNSMAX)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PDDC(0:KNSMAX,0:KNSMAX)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PDDD(0:KNSMAX,0:KNSMAX)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PDDE(0:KNSMAX,0:KNSMAX)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PDDH(0:KNSMAX)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PDDI(0:KNSMAX)

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZDLX,ZDLX1,ZDLSITA,ZDL1SITA,ZDLK,ZDL1,ZDLS,ZDLLDN

INTEGER(KIND=JPIM), PARAMETER :: JPKD=KIND(ZDLX)

INTEGER(KIND=JPIM) :: JM, JN, JK
REAL(KIND=JPRB) :: Z

!     ------------------------------------------------------------------

!*       1. First two columns.
!           ------------------

ZDLX=PDDMU
ZDLX1=ACOS(ZDLX)
ZDLSITA=SQRT(1.0_JPRB-ZDLX*ZDLX)

PDDPOL(0,0)=1._JPRB
ZDLLDN = 0.0_JPRB

! IF WE ARE LESS THAN 1Meter FROM THE POLE,
IF(ABS(REAL(ZDLSITA,KIND(Z))) <= SQRT(EPSILON(Z)))THEN
  ZDLX=1._JPRB
  ZDLSITA=0._JPRB
  ZDL1SITA=0._JPRB
ELSE
  ZDL1SITA=1.0_JPRB/ZDLSITA
ENDIF

!*          ordinary Legendre polynomials from series expansion
!           ---------------------------------------------------

! even N
DO JN=2,KNSMAX,2
  ZDLK = 0.5_JPRB*PFN(JN,0)
  ZDLLDN = 0.0_JPRB
  ! represented by only even k
  DO JK=2,JN,2
    ! normalised ordinary Legendre polynomial == \overbar{P_n}^0
    ZDLK = ZDLK + PFN(JN,JK)*COS(PDDI(JK)*ZDLX1)
    ! normalised associated Legendre polynomial == \overbar{P_n}^1
    ZDLLDN = ZDLLDN + PDDA(JN)*PFN(JN,JK)*PDDI(JK)*SIN(PDDI(JK)*ZDLX1)
  ENDDO
  PDDPOL(0,JN) = ZDLK
  PDDPOL(1,JN) = ZDLLDN
ENDDO
! odd N
DO JN=1,KNSMAX,2
  ZDLK = 0.0_JPRB
  ZDLLDN = 0.0_JPRB
  ! represented by only odd k
  DO JK=1,JN,2
    ! normalised ordinary Legendre polynomial == \overbar{P_n}^0
    ZDLK = ZDLK + PFN(JN,JK)*COS(PDDI(JK)*ZDLX1)
    ! normalised associated Legendre polynomial == \overbar{P_n}^1
    ZDLLDN = ZDLLDN + PDDA(JN)*PFN(JN,JK)*PDDI(JK)*SIN(PDDI(JK)*ZDLX1)
  ENDDO
  PDDPOL(0,JN) = ZDLK
  PDDPOL(1,JN) = ZDLLDN
ENDDO

!     ------------------------------------------------------------------

!*       2. Diagonal (the terms 0,0 and 1,1 have already been computed)
!           Belousov, equation (23)
!           -----------------------------------------------------------

ZDLS=ZDL1SITA*TINY(ZDLS)

#ifdef VPP
!OCL SCALAR
#endif
DO JN=2,KNSMAX
  PDDPOL(JN,JN)=PDDPOL(JN-1,JN-1)*ZDLSITA*PDDH(JN)
  IF ( ABS(PDDPOL(JN,JN))  <  ZDLS ) PDDPOL(JN,JN)=0.0_JPRB
ENDDO

!     ------------------------------------------------------------------

!*       3. General recurrence (Belousov, equation 17)
!           -----------------------------------------

DO JN=3,KNSMAX
!DIR$ IVDEP
!OCL NOVREC
  DO JM=2,JN-1
    PDDPOL(JM,JN)=PDDC(JM,JN)*PDDPOL(JM-2,JN-2)&
     &-PDDD(JM,JN)*PDDPOL(JM-2,JN-1)*ZDLX &
     &+PDDE(JM,JN)*PDDPOL(JM  ,JN-1)*ZDLX
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE SUPOL
END MODULE SUPOL_MOD
