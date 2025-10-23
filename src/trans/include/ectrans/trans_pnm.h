! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

INTERFACE
SUBROUTINE TRANS_PNM(KRESOL,KM,PRPNM,LDTRANSPOSE,LDCHEAP)

! begin_doc_block
! ## `TRANS_PNM`
!
! ### Signature
!
! ```f90
! SUBROUTINE TRANS_PNM(KRESOL, KM, PRPNM, LDTRANSPOSE, LDCHEAP)
! ```
!
! ### Purpose
!
! This subroutine computes the Legendre polynomials for a given wavenumber.
!
! ### `INTENT(IN)` arguments
!
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KM`  
!   The zonal wavenumber to compute the polynomials for.
!
! ### `OPTIONAL, INTENT(IN)` arguments
!
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KRESOL`  
!   Resolution handle returned by original call to `SETUP_TRANS`.  
!   *Default*: `1` (i.e. first resolution handle)
! - `LOGICAL, OPTIONAL, INTENT(IN) :: LDTRANSPOSE`  
!   Whether to transpose the output array (see dimensions of `PRPNM` below).  
!   *Default*: `.FALSE.`
! - `LOGICAL, OPTIONAL, INTENT(IN) :: LDCHEAP`  
!   Whether to use "cheap" algorithm to compute the polynomials (which are albeit less accurate).  
!   *Default*: `.FALSE.`
!
! ### `INTENT(OUT)` arguments
!
! - `REAL(KIND=JPRB), INTENT(OUT) :: PRPNM`  
!   The computed Legendre polynomials. The dimensions depend on the value of `LDTRANSPOSE`.  
!   Note that this must be already allocated.  
!   See [`SETUP_DIMS_MOD`](https://sites.ecmwf.int/docs/ectrans/sourcefile/setup_dims_mod.f90.html)
!   for the defininition of `R%NLEI3`.  
!   Dimensions: (`R%NTMAX-KM+3`, `R%NLEI3`) if `LDTRANSPOSE == .TRUE.`, else (`R%NLEI3`, `R%NTMAX-KM+3`)
! end_doc_block

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KRESOL
INTEGER(KIND=JPIM) ,INTENT(IN)  :: KM
REAL(KIND=JPRB)    ,INTENT(OUT) :: PRPNM(:,:)
LOGICAL, OPTIONAL, INTENT(IN) :: LDTRANSPOSE
LOGICAL, OPTIONAL, INTENT(IN) :: LDCHEAP

END SUBROUTINE TRANS_PNM
END INTERFACE
