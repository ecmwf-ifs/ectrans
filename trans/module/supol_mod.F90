MODULE SUPOL_MOD
CONTAINS
SUBROUTINE SUPOL(KNSMAX,PDDMU,PDDPOL, &
 & PDDA,PDDB,PDDC,PDDD,PDDE,PDDF,PDDG,PDDH,PDDI)

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
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE PARKIND2  ,ONLY : JPRH
#if defined(NECSX) && defined(REALHUGE)
USE QUAD_EMU_MOD
#endif

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KNSMAX
REAL(KIND=JPRH)   ,INTENT(IN)  :: PDDMU
REAL(KIND=JPRH)   ,INTENT(OUT) :: PDDPOL(0:KNSMAX,0:KNSMAX)
REAL(KIND=JPRH)   ,INTENT(IN)  :: PDDA(0:KNSMAX)
REAL(KIND=JPRH)   ,INTENT(IN)  :: PDDB(0:KNSMAX)
REAL(KIND=JPRH)   ,INTENT(IN)  :: PDDC(0:KNSMAX,0:KNSMAX)
REAL(KIND=JPRH)   ,INTENT(IN)  :: PDDD(0:KNSMAX,0:KNSMAX)
REAL(KIND=JPRH)   ,INTENT(IN)  :: PDDE(0:KNSMAX,0:KNSMAX)
REAL(KIND=JPRH)   ,INTENT(IN)  :: PDDF(0:KNSMAX)
REAL(KIND=JPRH)   ,INTENT(IN)  :: PDDG(0:KNSMAX)
REAL(KIND=JPRH)   ,INTENT(IN)  :: PDDH(0:KNSMAX)
REAL(KIND=JPRH)   ,INTENT(IN)  :: PDDI(0:KNSMAX)

!     ------------------------------------------------------------------

REAL(KIND=JPRH) :: ZDLX,ZDLSITA,ZDL1SITA,ZDLKM2,ZDLKM1,ZDLK,ZDL1,ZDLS

INTEGER(KIND=JPIM) :: JM, JN
REAL(KIND=JPRB) :: Z

#if defined(NECSX) && defined(REALHUGE)
REAL(KIND=JPRH), ALLOCATABLE, DIMENSION(:) :: ZDDPOL_
REAL(KIND=JPRH) :: ZDLXM
#endif

!     ------------------------------------------------------------------

!*       1. First two columns.
!           ------------------

ZDLX=PDDMU
ZDLSITA=SQRT(1.0_JPRB-ZDLX*ZDLX)

! IF WE ARE LESS THAN 1Meter FROM THE POLE,
IF(ABS(REAL(ZDLSITA,KIND(Z))) <= SQRT(EPSILON(Z)))THEN
  ZDLX=1._JPRB
  ZDLSITA=0._JPRB
  ZDL1SITA=0._JPRB
ELSE
  ZDL1SITA=1.0_JPRB/ZDLSITA
ENDIF
ZDLKM2=1._JPRB
ZDLKM1=ZDLX
PDDPOL(0,0)=ZDLKM2
PDDPOL(0,1)=ZDLKM1*PDDA(1)
PDDPOL(1,1)=ZDLSITA*PDDB(1)
DO JN=2,KNSMAX
  ZDLK=PDDF(JN)*ZDLX*ZDLKM1-PDDG(JN)*ZDLKM2
  ZDL1=PDDI(JN)*(ZDLKM1-ZDLX*ZDLK)*ZDL1SITA
  PDDPOL(0,JN)=ZDLK*PDDA(JN)
  PDDPOL(1,JN)=ZDL1*PDDB(JN)
  ZDLKM2=ZDLKM1
  ZDLKM1=ZDLK
ENDDO

!     ------------------------------------------------------------------

!*       2. Diagonal (the terms 0,0 and 1,1 have already been computed)
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

!*       3. General recurrence.
!           -------------------

#if defined(NECSX) && defined(REALHUGE)

ALLOCATE(ZDDPOL_(2:KNSMAX-1),STAT=JN)

ZDLXM = -ZDLX
DO JN=3,KNSMAX
  CALL QXMY(PDDC(2:JN-1,JN),PDDPOL(0:JN-3,JN-2),JN-2,PDDPOL (2:JN-1,JN))
  CALL QXMY(PDDD(2:JN-1,JN),PDDPOL(0:JN-3,JN-1),JN-2,ZDDPOL_(2:JN-1))
  CALL QAXPY(ZDLXM,ZDDPOL_(2),PDDPOL(2,JN),JN-2)
  CALL QXMY(PDDE(2:JN-1,JN),PDDPOL(2:JN-1,JN-1),JN-2,ZDDPOL_(2:JN-1))
  CALL QAXPY(ZDLX,ZDDPOL_(2),PDDPOL(2,JN),JN-2)
ENDDO

DEALLOCATE (ZDDPOL_)

#else

DO JN=3,KNSMAX
!DIR$ IVDEP
!OCL NOVREC
  DO JM=2,JN-1
    PDDPOL(JM,JN)=PDDC(JM,JN)*PDDPOL(JM-2,JN-2)&
     & -PDDD(JM,JN)*PDDPOL(JM-2,JN-1)*ZDLX &
     & +PDDE(JM,JN)*PDDPOL(JM  ,JN-1)*ZDLX
  ENDDO
ENDDO

#endif

!     ------------------------------------------------------------------

END SUBROUTINE SUPOL
END MODULE SUPOL_MOD
