! (C) Copyright 2025- Meteo-France.
! (C) Copyright 2025- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

SUBROUTINE GET_LEGENDRE_ASSETS(KRETURNCODE, NLAT, KTRUNC, KSPOLEGL, KLOEN, KNMENG, PGW, PRPNM)
! ** PURPOSE
!    Simplified wrapper to TRANS_INQ for obtaining arrays necessary for performing Legendre
!    transform (Gaussian weights, Legendre polynomials and NMENG (cutoff zonal wavenumber for each
!    latitude))
!
! ** DUMMY ARGUMENTS
!    NLAT: number of latitudes (pole-to-pole) in grid-point space
!    KTRUNC: truncation
!    KSPOLEGL: Size of second dimension of Legendre polynomials
!    KLOEN: number of points on each latitude row
!    KNMENG: cut-off zonal wavenumber
!    PGW: Gaussian weights
!    PRPNM: Legendre polynomials
!
! ** AUTHOR
!    2 July 2025, S. Hatfield

! I. Dummy arguments declaration
USE ISO_FORTRAN_ENV, ONLY: INT32, REAL64
IMPLICIT NONE
INTEGER(KIND=INT32), INTENT(OUT) :: KRETURNCODE
INTEGER(KIND=INT32), INTENT(IN) :: NLAT
INTEGER(KIND=INT32), INTENT(IN) :: KTRUNC
INTEGER(KIND=INT32), INTENT(IN) :: KSPOLEGL
INTEGER(KIND=INT32), DIMENSION(NLAT), INTENT(IN) :: KLOEN
INTEGER(KIND=INT32), DIMENSION(NLAT), INTENT(OUT) :: KNMENG
REAL(KIND=REAL64), DIMENSION(NLAT), INTENT(OUT) :: PGW
REAL(KIND=REAL64), DIMENSION(NLAT/2,KSPOLEGL), INTENT(OUT) :: PRPNM

! II. Local variables declaration
INTEGER, DIMENSION(NLAT) :: INMENG
REAL(KIND=REAL64), DIMENSION(NLAT) :: ZGW
REAL(KIND=REAL64), DIMENSION(NLAT/2,KSPOLEGL) :: ZRPNM

#include "setup_trans0.h"
#include "setup_trans.h"
#include "trans_inq.h"
#include "trans_end.h"

CALL SETUP_TRANS0(LDMPOFF=.TRUE.)
CALL SETUP_TRANS(KSMAX=KTRUNC, KDGL=NLAT, KLOEN=KLOEN)
CALL TRANS_INQ(KNMENG=INMENG, PGW=ZGW, PRPNM=ZRPNM)
KNMENG = INMENG
PGW = ZGW
PRPNM = ZRPNM
CALL TRANS_END

END SUBROUTINE GET_LEGENDRE_ASSETS
