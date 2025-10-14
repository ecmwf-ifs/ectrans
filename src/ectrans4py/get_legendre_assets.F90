! (C) Copyright 2025- Meteo-France.
! (C) Copyright 2025- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

SUBROUTINE GET_LEGENDRE_ASSETS(KRETURNCODE, NLAT, KTRUNC, KSPOLEGL, KLOEN, KNUMMAXRESOL, &
  &                            KNMENG, PGW, PRPNM)
! ** PURPOSE
!    Simplified wrapper to TRANS_INQ for obtaining arrays necessary for performing Legendre transform
!    (Gaussian weights, Legendre polynomials and NMENG (cutoff zonal wavenumber for each latitude))
!
! ** DUMMY ARGUMENTS
!    NLAT: number of latitudes (pole-to-pole) in grid-point space
!    KTRUNC: truncation
!    KSPOLEGL: Size of second dimension of Legendre polynomials
!    KLOEN: number of points on each latitude row
!    KNUMMAXRESOL: maximum number of troncatures handled
!    KNMENG: cut-off zonal wavenumber
!    PGW: Gaussian weights
!    PRPNM: Legendre polynomials
!
! ** AUTHOR
!    2 July 2025, S. Hatfield

! I. Dummy arguments declaration
USE ISO_FORTRAN_ENV, ONLY: INT64, REAL64
IMPLICIT NONE
INTEGER(KIND=INT64), INTENT(OUT) :: KRETURNCODE
INTEGER(KIND=INT64), INTENT(IN) :: NLAT
INTEGER(KIND=INT64), INTENT(IN) :: KTRUNC
INTEGER(KIND=INT64), INTENT(IN) :: KSPOLEGL
INTEGER(KIND=INT64), DIMENSION(NLAT), INTENT(IN) :: KLOEN
INTEGER(KIND=INT64), INTENT(IN) :: KNUMMAXRESOL
INTEGER(KIND=INT64), DIMENSION(NLAT), INTENT(OUT) :: KNMENG
REAL(KIND=REAL64), DIMENSION(NLAT), INTENT(OUT) :: PGW
REAL(KIND=REAL64), DIMENSION(NLAT/2,KSPOLEGL), INTENT(OUT) :: PRPNM

! II. Local variables declaration
INTEGER, DIMENSION(NLAT) :: ILOEN
INTEGER :: ISIZEI, ISIZEJ, IPHYSICALSIZEI, IPHYSICALSIZEJ, ITRUNCX, ITRUNCY, INUMMAXRESOL
LOGICAL :: LLSTOP
INTEGER :: IIDENTRESOL
INTEGER, DIMENSION(NLAT) :: INMENG
REAL(KIND=REAL64), DIMENSION(NLAT) :: ZGW
REAL(KIND=REAL64), DIMENSION(NLAT/2,KSPOLEGL) :: ZRPNM
REAL(KIND=REAL64) :: ZDELTAX, ZDELTAY

#include "trans_inq.h"

ILOEN(:) = KLOEN(:)
ISIZEI = 0
ISIZEJ = NLAT
IPHYSICALSIZEI = 0
IPHYSICALSIZEJ = 0
ITRUNCX = KTRUNC
ITRUNCY = 0
INUMMAXRESOL = KNUMMAXRESOL

! III. Setup
ZDELTAX = 0.0_REAL64
ZDELTAY = 0.0_REAL64
CALL SPEC_SETUP4PY(KRETURNCODE, ISIZEI, ISIZEJ, IPHYSICALSIZEI, IPHYSICALSIZEJ, &
                  &ITRUNCX, ITRUNCY, INUMMAXRESOL, ILOEN, .FALSE., SIZE(ILOEN), &
                  &ZDELTAX, ZDELTAY, IIDENTRESOL, LLSTOP)
IF (.NOT. LLSTOP) THEN
  CALL TRANS_INQ(KRESOL=IIDENTRESOL, KNMENG=INMENG, PGW=ZGW, PRPNM=ZRPNM)
  KNMENG = INMENG
  PGW = ZGW
  PRPNM = ZRPNM
ENDIF

END SUBROUTINE GET_LEGENDRE_ASSETS
