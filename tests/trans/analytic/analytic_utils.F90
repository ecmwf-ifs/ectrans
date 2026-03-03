MODULE ANALYTIC_SOLUTIONS_MOD

  USE PARKIND1, ONLY: JPIM, JPRD

  REAL(KIND=JPRD) :: Z_PI = 4.0_JPRD * ATAN(1.0_JPRD)
  REAL(KIND=JPRD), ALLOCATABLE :: ZMU(:), LEGPOLYS(:, :, :), GELAM(:, :), GELAT(:, :)
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
        NLATIDXS(JROF, IBL) = JGLAT
        JROF = JROF + 1
        IF (JROF > KPROMA) THEN
          JROF = 1
          IBL  = IBL + 1
        ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE ANALYTIC_INIT

  !===================================================================================================  
  ! Load the Legendre polynomials into the helper array LEGPOLYS
  !===================================================================================================

  SUBROUTINE BUFFER_LEGENDRE_POLYNOMIALS(NSMAX, NDGL)

    USE PARKIND1, ONLY: JPRB

    IMPLICIT NONE

    INTEGER(KIND=JPIM), INTENT(IN) :: NSMAX
    INTEGER(KIND=JPIM), INTENT(IN) :: NDGL

    REAL(KIND=JPRB), DIMENSION(NDGL / 2, (NSMAX + 2) * (NSMAX + 3) / 2) :: RPNM
    REAL(KIND=JPRB), DIMENSION(NDGL, (NSMAX + 2) * (NSMAX + 3) / 2) :: RPNM_FULL
    INTEGER(KIND=JPIM) :: JNM, JN, JNINV, JM, ILAT

    ALLOCATE(LEGPOLYS(NFIRSTLAT:NLASTLAT, 0:NSMAX, 0:NSMAX))
    LEGPOLYS = 0.0

    CALL TRANS_INQ(PRPNM=RPNM)

    DO ILAT = 1, NDGL / 2
      JNM = 0
      DO JM = 0, NSMAX + 1
        DO JNINV = JM, NSMAX + 1
          JN = NSMAX + 1 - JNINV + JM
          JNM = JNM + 1
          IF (MOD(JM + JN, 2) == 0) THEN
            RPNM_FULL(ILAT,            JNM) = RPNM(ILAT, JNM)
            RPNM_FULL(NDGL - ILAT + 1, JNM) = RPNM(ILAT, JNM)
          ELSE
            RPNM_FULL(ILAT,            JNM) = RPNM(ILAT, JNM)
            RPNM_FULL(NDGL - ILAT + 1, JNM) = -RPNM(ILAT, JNM)
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    DO ILAT = NFIRSTLAT, NLASTLAT
      JNM = 0
      DO JM = 0, NSMAX + 1
        DO JNINV = JM, NSMAX + 1
          JN = NSMAX + 1 - JNINV + JM
          JNM = JNM + 1
          IF (JM <= NSMAX .AND. JN <= NSMAX) THEN
            LEGPOLYS(ILAT, JM, JN) = RPNM_FULL(ILAT - NFIRSTLAT + 1, JNM)
          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE BUFFER_LEGENDRE_POLYNOMIALS

  !===================================================================================================  
  ! Deallocate the helper arrays used for the analytic solutions.
  !===================================================================================================

  SUBROUTINE ANALYTIC_END()

    IMPLICIT NONE

    DEALLOCATE(ZMU, NLATIDXS, GELAM, GELAT)
    DEALLOCATE(LEGPOLYS)

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
          PSPH_ANALYTIC(JROF, IBL) = LEGPOLYS(NLATIDXS(JROF, IBL), KZONAL, KTOTAL)
        ELSE
          PSPH_ANALYTIC(JROF, IBL) = 2.0_JPRD * COS(KZONAL * GELAM(JROF, IBL)) * &
            & LEGPOLYS(NLATIDXS(JROF, IBL), KZONAL, KTOTAL)
        ENDIF
      ENDDO
    ENDDO

  END FUNCTION COMPUTE_ANALYTIC_SOLUTION

  !===================================================================================================

END MODULE ANALYTIC_SOLUTIONS_MOD