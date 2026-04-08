MODULE ANALYTIC_UTILS

  USE PARKIND1, ONLY: JPIM, JPRD, JPRB

  REAL(KIND=JPRD) :: Z_PI = 4.0_JPRD * ATAN(1.0_JPRD)
  REAL(KIND=JPRD), ALLOCATABLE :: ZMU(:), GELAM(:, :), GELAT(:, :)
  REAL(KIND=JPRD), ALLOCATABLE :: PLEGPOLYS(:, :)
  INTEGER(KIND=JPIM), ALLOCATABLE :: NLATIDXS(:, :)
  INTEGER(KIND=JPIM) :: NFIRSTLAT, NLASTLAT

#include "trans_inq.h"

  CONTAINS

  !===================================================================================================
  ! Compute with the help of TRANS_INQ the geographic longitude GELAM and latitude GELAM.
  ! Also create a helper array NLATIDXS(NPROMA, NGPBLKS) which contains for each blocked point the
  ! global latitude index. This is used later to retrieve the corresponding Legendre polynomial.
  !===================================================================================================

  SUBROUTINE ANALYTIC_INIT(KPROMA, KGPBLKS, KDGL, K_REGIONS_NS, K_REGIONS_EW, KLOEN)

    IMPLICIT NONE

    INTEGER(KIND=JPIM), INTENT(IN) :: KPROMA
    INTEGER(KIND=JPIM), INTENT(IN) :: KGPBLKS
    INTEGER(KIND=JPIM), INTENT(IN) :: KDGL
    INTEGER(KIND=JPIM), INTENT(IN) :: K_REGIONS_NS
    INTEGER(KIND=JPIM), INTENT(IN) :: K_REGIONS_EW
    INTEGER(KIND=JPIM), DIMENSION(KDGL), INTENT(IN) :: KLOEN

    INTEGER(KIND=JPIM) :: NPTRFLOFF, MY_REGION_NS, MY_REGION_EW
    INTEGER(KIND=JPIM) :: JGLAT, IOFF, ILAT, ISTLON, IENDLON, JLON, JROF, IBL
    INTEGER(KIND=JPIM), DIMENSION(K_REGIONS_NS) :: NFRSTLAT, NLSTLAT
    INTEGER(KIND=JPIM), DIMENSION(KDGL + K_REGIONS_NS - 1, K_REGIONS_EW) :: NSTA, NONL
    REAL(KIND=JPRD) :: ZLAT, ZLON

    ALLOCATE(ZMU(KDGL), NLATIDXS(KPROMA, KGPBLKS), GELAM(KPROMA, KGPBLKS), GELAT(KPROMA, KGPBLKS))

    CALL TRANS_INQ(KPTRFLOFF=NPTRFLOFF, KMY_REGION_NS=MY_REGION_NS, KMY_REGION_EW=MY_REGION_EW, KFRSTLAT=NFRSTLAT, &
      &            KLSTLAT=NLSTLAT, KSTA=NSTA, KONL=NONL, PMU=ZMU)

    ILAT = NPTRFLOFF
    IBL  = 1
    JROF = 1
    NFIRSTLAT = NFRSTLAT(MY_REGION_NS)
    NLASTLAT = NLSTLAT(MY_REGION_NS)
    DO JGLAT = NFIRSTLAT, NLASTLAT
      ZLAT = ASIN(ZMU(JGLAT))
      ILAT = ILAT + 1
      ISTLON = NSTA(ILAT, MY_REGION_EW)
      IENDLON = ISTLON - 1 + NONL(ILAT, MY_REGION_EW)
      DO JLON = ISTLON, IENDLON
        ZLON = REAL(JLON - 1, JPRD) * 2.0_JPRD * Z_PI / REAL(KLOEN(JGLAT), JPRD)
        GELAM(JROF, IBL) = ZLON ! Longitude
        GELAT(JROF, IBL) = ZLAT ! Latitude
        NLATIDXS(JROF, IBL) = JGLAT - NFIRSTLAT + 1 ! Latitude of this (JROF, IBL) relative to NFIRSTLAT
        JROF = JROF + 1
        IF (JROF > KPROMA) THEN
          JROF = 1
          IBL  = IBL + 1
        ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE ANALYTIC_INIT

  !===================================================================================================
  ! Compute the Legendre polynomials for all latitudes and total wavenumbers at the given zonal
  ! wavenumber.
  !===================================================================================================

  SUBROUTINE PREPARE_LEGENDRE_POLYNOMIALS(KZONAL, KSMAX)

    USE TPM_POL, ONLY: INI_POL, END_POL
    USE SUPOLF_MOD, ONLY: SUPOLF

    IMPLICIT NONE

    INTEGER(KIND=JPIM), INTENT(IN) :: KZONAL
    INTEGER(KIND=JPIM), INTENT(IN) :: KSMAX

    INTEGER(KIND=JPIM) :: JGLAT

    IF (ALLOCATED(PLEGPOLYS)) DEALLOCATE(PLEGPOLYS)
    ALLOCATE(PLEGPOLYS(NLASTLAT - NFIRSTLAT + 1, 0:KSMAX))

    CALL INI_POL(KSMAX)

    DO JGLAT = NFIRSTLAT, NLASTLAT
      CALL SUPOLF(KZONAL, KSMAX, ZMU(JGLAT), PLEGPOLYS(JGLAT - NFIRSTLAT + 1, :))
    ENDDO

    CALL END_POL

  END SUBROUTINE PREPARE_LEGENDRE_POLYNOMIALS

  !===================================================================================================
  ! Deallocate the helper arrays used for the analytic solutions.
  !===================================================================================================

  SUBROUTINE ANALYTIC_END()

    IMPLICIT NONE

    DEALLOCATE(ZMU, NLATIDXS, GELAM, GELAT)
    DEALLOCATE(PLEGPOLYS)

  END SUBROUTINE ANALYTIC_END

  !====================================================================================================
  ! Compute analytic solution for a specific total wavenumber n and zonal wavenumber m by going through
  ! all points and using the point-wise function analytic_spherical_harmonic_point.
  !====================================================================================================

  FUNCTION COMPUTE_ANALYTIC_SOLUTION(KPROMA, KGPBLKS, KGPTOT, KZONAL, KTOTAL) RESULT(PSPH_ANALYTIC)

    IMPLICIT NONE

    INTEGER(KIND=JPIM), INTENT(IN) :: KPROMA
    INTEGER(KIND=JPIM), INTENT(IN) :: KGPBLKS
    INTEGER(KIND=JPIM), INTENT(IN) :: KGPTOT
    INTEGER(KIND=JPIM), INTENT(IN) :: KZONAL
    INTEGER(KIND=JPIM), INTENT(IN) :: KTOTAL
    REAL(KIND=JPRD), DIMENSION(KPROMA, KGPBLKS) :: PSPH_ANALYTIC

    INTEGER(KIND=JPIM) :: JKGLO, IEND, IBL, JROF

    DO JKGLO = 1, KGPTOT, KPROMA
      IEND = MIN(KPROMA, KGPTOT - JKGLO + 1)
      IBL  = (JKGLO - 1) / KPROMA + 1
      DO JROF = 1, IEND
        IF (KZONAL == 0) THEN
          PSPH_ANALYTIC(JROF, IBL) = PLEGPOLYS(NLATIDXS(JROF, IBL), KTOTAL)
        ELSE
          PSPH_ANALYTIC(JROF, IBL) = 2.0_JPRD * COS(KZONAL * GELAM(JROF, IBL)) * &
            & PLEGPOLYS(NLATIDXS(JROF, IBL), KTOTAL)
        ENDIF
      ENDDO
    ENDDO

  END FUNCTION COMPUTE_ANALYTIC_SOLUTION

  !====================================================================================================
  ! Compute the maximum error between two blocked grid point arrays.
  ! This is needed when comparing arrays that have a partially filled <NPROMA tail block.
  !====================================================================================================

  FUNCTION COMPUTE_MAX_ERROR(KGPTOT, KPROMA, PGP1, PGP2) RESULT(MAX_ERROR)

    IMPLICIT NONE

    INTEGER(KIND=JPIM), INTENT(IN) :: KGPTOT
    INTEGER(KIND=JPIM), INTENT(IN) :: KPROMA
    REAL(KIND=JPRB), DIMENSION(:, :), INTENT(IN) :: PGP1
    REAL(KIND=JPRD), DIMENSION(:, :), INTENT(IN) :: PGP2
    REAL(KIND=JPRD) :: MAX_ERROR

    INTEGER(KIND=JPIM) :: JKGLO, IEND, IBL, JROF

    MAX_ERROR = 0.0_JPRD

    DO JKGLO = 1, KGPTOT, KPROMA
      IEND = MIN(KPROMA, KGPTOT - JKGLO + 1)
      IBL = (JKGLO - 1) / KPROMA + 1
      DO JROF = 1, IEND
        MAX_ERROR = MAX(MAX_ERROR, ABS(PGP1(JROF, IBL) - PGP2(JROF, IBL)))
      ENDDO
    ENDDO
  END FUNCTION COMPUTE_MAX_ERROR

  !===================================================================================================

END MODULE ANALYTIC_UTILS