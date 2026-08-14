! (C) Copyright 1992- ECMWF.
! (C) Copyright 1992- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE GAWL_MOD
CONTAINS
SUBROUTINE GAWL(PFN, PL, PW, PEPS, KN, KITER, PMOD)

!**** *GAWL * - Routine to perform the Newton loop

!     Purpose.
!     --------
!           Find 0 of Legendre polynomial with Newton loop
!**   Interface.
!     ----------
!        *CALL* *GAWL(PFN,PL,PW,PEPS,KN,KITER,PMOD)

!        Explicit arguments :
!        --------------------
! PFN    Fourier coefficients of series expansion
!        for the ordinary Legendre polynomials     (in)
! PL     Gaussian latitude                         (inout)
! PW     Gaussian weight                           (out)
! PEPS   0 of the machine                          (in)
! KN     Truncation                                (in)
! KITER  Number of iterations                      (out)
! PMOD   Last modification                         (inout)

!        Implicit arguments :
!        --------------------
!       None

!     Method.
!     -------
!        Newton Loop.

!     Externals.
!     ----------
!        CPLEDN

!     Reference.
!     ----------

!     ARPEGE Documentation vol.2, ch3.

!     Author.
!     -------
!        Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 92-12-18
!        K. Yessad (Sep 2008): cleaning, improve comments.
!        Nils Wedi + Mats Hamrud, 2009-02-05 revised following Swarztrauber, 2002
!      F. Vana  05-Mar-2015  Support for single precision
!     ------------------------------------------------------------------

USE EC_PARKIND, ONLY: JPRD, JPIM

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)    :: KN
REAL(KIND=JPRD),    INTENT(IN)    :: PFN(0:KN/2)
REAL(KIND=JPRD),    INTENT(INOUT) :: PL
REAL(KIND=JPRD),    INTENT(OUT)   :: PW
REAL(KIND=JPRD),    INTENT(IN)    :: PEPS
INTEGER(KIND=JPIM), INTENT(OUT)   :: KITER
REAL(KIND=JPRD),    INTENT(INOUT) :: PMOD

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: ITEMAX = 20
INTEGER(KIND=JPIM) :: IODD
REAL(KIND=JPRD) :: ZDLK, ZDLLDN
INTEGER(KIND=JPIM) :: IK, JN

!     ------------------------------------------------------------------

!*       1. Initialization.
!           ---------------

IODD = MOD(KN, 2)

!     ------------------------------------------------------------------

!*       2. Newton iteration.
!           -----------------

PMOD = HUGE(1.0_JPRD)

DO KITER = 1, ITEMAX + 1
  ZDLLDN = 0.0_JPRD
  IK = 1
  IF (ABS(PMOD) <= PEPS * 1000._JPRD) THEN
    ! Last pass
    DO JN = 2 - IODD, KN, 2
      ! normalised derivative
      ZDLLDN = ZDLLDN - PFN(IK) * REAL(JN, JPRD) * SIN(REAL(JN, JPRD) * PL)
      IK = IK + 1
    ENDDO
    PW = REAL(2 * KN + 1, JPRD) / ZDLLDN ** 2
    EXIT
  ENDIF

  ZDLK = 0.0_JPRD
  IF (IODD == 0) ZDLK = 0.5_JPRD * PFN(0)

  DO JN = 2 - IODD, KN, 2
    ! normalised ordinary Legendre polynomial == \overbar{P_n}^0
    ZDLK = ZDLK + PFN(IK) * COS(REAL(JN, JPRD) * PL)
    ! normalised derivative == d/d\theta(\overbar{P_n}^0)
    ZDLLDN = ZDLLDN - PFN(IK) * REAL(JN, JPRD) * SIN(REAL(JN, JPRD) * PL)
    IK = IK + 1
  ENDDO
  ! Newton method
  PMOD = -ZDLK / ZDLLDN
  PL = PL + PMOD
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE GAWL
END MODULE GAWL_MOD
