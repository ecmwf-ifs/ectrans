! (C) Copyright 1987- ECMWF.
! (C) Copyright 1987- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE SUPOLF_MOD
CONTAINS
SUBROUTINE SUPOLF(KM,KNSMAX,DDMU,DDPOL,KCHEAP)

!**** *SUPOL * - Routine to compute the Legendre polynomials

!     Purpose.
!     --------
!           For a given value of mu and M, computes the Legendre
!           polynomials upto KNSMAX

!**   Interface.
!     ----------
!        *CALL* *SUPOLF(KM,KNSMAX,DDMU,DDPOL,KCHEAP)

!        Explicit arguments :
!        --------------------
!              KM       :  zonal wavenumber M
!              KNSMAX   :  Truncation  (triangular)
!              DDMU     :  Abscissa at which the polynomials are computed (mu)
!              DDPOL    :  Polynomials (the first index is m and the second n)
!              KCHEAP   :  odd/even saving switch


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
!        Nils Wedi + George Mozdzynski + Mats Hamrud

!     Modifications.
!     --------------
!        Original : 87-10-15
!        K. YESSAD (MAY 1998): modification to avoid underflow.
!        R. El Khatib 11-Apr-2007 Emulation of vectorized quadruple precision
!                                 on NEC
!      F. Vana  05-Mar-2015  Support for single precision
!     ------------------------------------------------------------------

USE EC_PARKIND  ,ONLY : JPRD, JPIM

USE TPM_POL   ,ONLY : DFI, DFB, DFG, DFA, DFF

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KM
INTEGER(KIND=JPIM),INTENT(IN)  :: KNSMAX
REAL(KIND=JPRD)   ,INTENT(IN)  :: DDMU
REAL(KIND=JPRD)   ,INTENT(OUT) :: DDPOL(0:KNSMAX)

INTEGER(KIND=JPIM),INTENT(IN), OPTIONAL :: KCHEAP

REAL(KIND=JPRD) :: DLX, ZCOS_THETA, ZCOS2_THETA, ZCOS_THETA_R, DLK, DL1, DLKM1, DLKM2

INTEGER(KIND=JPIM) :: JN, KKL, ICHEAP, IC, IEND
REAL(KIND=JPRD) :: DCL, DDL

REAL(KIND=JPRD) :: ZFAC, ZLSITA, ZFAC0, ZFAC1, ZMULT, ZEPS

INTEGER(KIND=JPIM) :: JCORR, ICORR3, ICORR(KNSMAX)
REAL(KIND=JPRD) :: ZSCALE, ZISCALE

DCL(KKL)=SQRT((REAL(KKL-KM+1,JPRD)*REAL(KKL-KM+2,JPRD)* &
 & REAL(KKL+KM+1,JPRD)*REAL(KKL+KM+2,JPRD))/(REAL(2*KKL+1,JPRD)*REAL(2*KKL+3,JPRD)*&
 & REAL(2*KKL+3,JPRD)*REAL(2*KKL+5,JPRD)))
DDL(KKL)=(2.0_JPRD*REAL(KKL,JPRD)*REAL(KKL+1,JPRD)-2.0_JPRD*REAL(KM**2,JPRD)-1.0_JPRD)/ &
 & (REAL(2*KKL-1,JPRD)*REAL(2*KKL+3,JPRD))

!     ------------------------------------------------------------------

!*       1. First two columns.
!           ------------------

ZEPS  = EPSILON(ZSCALE)
ICORR3 = 0

ICHEAP = 1
IF (PRESENT(KCHEAP)) THEN
  ICHEAP = KCHEAP
ENDIF

! Ordinate to compute polynomials at (i.e. the sine of latitude)
DLX = DDMU

! Cosine^2 of latitude
ZCOS2_THETA = 1.0_JPRD - DLX * DLX
ZCOS_THETA = SQRT(ZCOS2_THETA)

!*          ordinary Legendre polynomials from series expansion
!           ---------------------------------------------------

! this is supol_fast just using single KM
IF (ABS(ZCOS_THETA) <= ZEPS) THEN
  DLX = 1.0_JPRD
  ZCOS_THETA = 0.0_JPRD
  ZCOS_THETA_R = 0.0_JPRD
  ZCOS2_THETA = 0.0_JPRD
ELSE
  ZCOS_THETA_R = 1.0_JPRD / ZCOS_THETA
ENDIF

DLKM2 = 1._JPRD
DLKM1 = DLX

IF (KM == 0) THEN
  DDPOL(0) = DLKM2
  DDPOL(1) = DLKM1 * DFB(1) / DFA(1)
  DO JN = 2, KNSMAX
    DLK = DFF(JN) * DLX * DLKM1 - DFG(JN) * DLKM2
    DDPOL(JN) = DLK * DFB(JN) / DFA(JN)
    DLKM2 = DLKM1
    DLKM1 = DLK
  ENDDO
ELSEIF (KM == 1) THEN
  DDPOL(0) = 0
  DDPOL(1) = ZCOS_THETA * DFB(1)
  DO JN = 2, KNSMAX
    DLK = DFF(JN) * DLX * DLKM1 - DFG(JN) * DLKM2
    DL1 = DFI(JN) * (DLKM1 - DLX * DLK) * ZCOS_THETA_R
    DDPOL(JN) = DL1 * DFB(JN)
    DLKM2 = DLKM1
    DLKM1 = DLK
  ENDDO
ELSE

!     ------------------------------------------------------------------
!*       KM >= 2
!     ------------------------------------------------------------------

  ! Maintaining the consistency with the CY41R1 reference
  ZSCALE=1.0E+100_JPRD
  ZISCALE=1.0E-100_JPRD

  ! Calculate (cos^2(theta))^(m/2) = (1 - mu^2)^(m/2)
  IEND = KM / 2
  ZLSITA = 1.0_JPRD
  DO JN = 1, IEND
    ZLSITA = ZLSITA * ZCOS2_THETA
    IF (ABS(ZLSITA) < ZISCALE) THEN
      ZLSITA = ZLSITA * ZSCALE
      ICORR3 = ICORR3 + 1
    ENDIF
  ENDDO
  IF (MOD(KM,2) == 1) ZLSITA = ZLSITA * ZCOS_THETA

  ! Calculate the first factorial term p_0 = sqrt(2m-1) * prod_{n=1}^{m-1} sqrt((2n-1)/(2n))
  ZFAC = 1.0_JPRD
  DO JN = 1, KM - 1
    ZFAC = ZFAC * SQRT(REAL(2 * JN - 1, JPRD))
    ZFAC = ZFAC / SQRT(REAL(2 * JN, JPRD))
  ENDDO
  ZFAC = ZFAC * SQRT(REAL(2 * KM - 1, JPRD))

  ZFAC0 = 1.0_JPRD
  ! Fill the first 4 values using explicit formulae
  DO IC = 0, MIN(KNSMAX - KM, 3)
    ! (2m+i)!
    ZFAC0 = ZFAC0 * REAL(2 * KM + IC, JPRD)

    SELECT CASE (IC)
    CASE (0)
      ZFAC1 = 1.0_JPRD ! 0!
      ! d_0 = d_0
      ZMULT = ZFAC
    CASE (1)
      ZFAC1 = 1.0_JPRD ! 1!
      ! p_1 = (2m+1) * p_0
      ZFAC = ZFAC * REAL(2 * KM + 1, JPRD)
      ! d_1 = mu * p_1
      ZMULT = ZFAC * DLX
    CASE (2)
      ZFAC1 = 2.0_JPRD ! 2!
      ! d_2 = 0.5 * p_2 * ( (2m+3) * mu^2 - 1 ) (p_2 = p_1)
      ZMULT = 0.5_JPRD * ZFAC * (REAL(2 * KM + 3, JPRD) * DLX * DLX - 1.0_JPRD)
    CASE (3)
      ZFAC1 = 6.0_JPRD ! 3!
      ! p_3 = (2m+3) * p_2
      ZFAC = ZFAC * REAL(2 * KM + 3, JPRD)
      ! d_3 = (1/6) * mu * p_3 * ( (2m+5) * mu^2 - 3 )
      ZMULT = (1.0_JPRD / 6.0_JPRD) * DLX * ZFAC * (REAL(2 * KM + 5, JPRD) * DLX * DLX - 3.0_JPRD)
    END SELECT

    ! P_{m,m+j} = (1 - mu^2)^{m/2} * d_j * sqrt( (2(m+j)+1) * j! / prod_{j'=0}^{j} (2m+j') )
    DDPOL(KM + IC) = ZLSITA * ZMULT * SQRT(2.0_JPRD * (REAL(KM + IC, JPRD) + 0.5_JPRD) * ZFAC1 / ZFAC0)
  ENDDO

  ICORR(:)=ICORR3
  IF (ICHEAP == 2) THEN
    ! symmetric case
    DO JN = KM + 2, KNSMAX - 2, 2
      IF (ABS(DDPOL(JN-2)) > ZSCALE) THEN
        DDPOL(JN-2) = DDPOL(JN-2) / ZSCALE
        DDPOL(JN) = DDPOL(JN) / ZSCALE
        ICORR(JN-2:KNSMAX) = ICORR(JN-2:KNSMAX) - 1
      ENDIF

      DDPOL(JN+2) = ((DLX * DLX - DDL(JN)) * DDPOL(JN) - DCL(JN-2) * DDPOL(JN-2)) / DCL(JN)
    ENDDO

    DO JN = KM, KNSMAX, 2
      DO JCORR = 1, ICORR(JN)
        DDPOL(JN) = DDPOL(JN) / ZSCALE
        IF (DDPOL(JN) < ZEPS) THEN
          DDPOL(JN) = ZEPS
        ENDIF
      ENDDO
    ENDDO

  ELSEIF (ICHEAP == 3) THEN
    ! antisymmetric case
    DO JN = KM + 3, KNSMAX - 2, 2
      IF (ABS(DDPOL(JN-2)) > ZSCALE) THEN
        DDPOL(JN-2) = DDPOL(JN-2) / ZSCALE
        DDPOL(JN) = DDPOL(JN) / ZSCALE
        ICORR(JN-2:KNSMAX) = ICORR(JN-2:KNSMAX) - 1
      ENDIF

      DDPOL(JN+2) = ((DLX * DLX - DDL(JN)) * DDPOL(JN) - DCL(JN-2) * DDPOL(JN-2)) / DCL(JN)
    ENDDO

    DO JN = KM + 1, KNSMAX, 2
      DO JCORR = 1, ICORR(JN)
        DDPOL(JN) = DDPOL(JN) / ZSCALE
        IF (DDPOL(JN) < ZEPS) THEN
          DDPOL(JN) = ZEPS
        ENDIF
      ENDDO
    ENDDO
  ELSE
    DO JN = KM + 2, KNSMAX - 2
      IF (ABS(DDPOL(JN-2)) > ZSCALE) THEN
        DDPOL(JN-2) = DDPOL(JN-2) / ZSCALE
        DDPOL(JN-1) = DDPOL(JN-1) / ZSCALE
        DDPOL(JN) = DDPOL(JN) / ZSCALE
        DDPOL(JN+1) = DDPOL(JN+1) / ZSCALE
        ICORR(JN-2:KNSMAX) = ICORR(JN-2:KNSMAX) - 1
      ENDIF

      DDPOL(JN+2) = ((DLX * DLX - DDL(JN)) * DDPOL(JN) - DCL(JN-2) * DDPOL(JN-2)) / DCL(JN)
    ENDDO

    DO JN = KM, KNSMAX
      DO JCORR = 1, ICORR(JN)
        DDPOL(JN) = DDPOL(JN) / ZSCALE
        IF (DDPOL(JN) < ZEPS) THEN
          DDPOL(JN) = ZEPS
        ENDIF
      ENDDO
    ENDDO
  ENDIF
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE SUPOLF
END MODULE SUPOLF_MOD
