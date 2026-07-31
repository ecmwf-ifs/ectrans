! (C) Copyright 2024- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

INTERFACE
SUBROUTINE GPNORM_TRANSTL(PGP,KFIELDS,KPROMA,PAVE,KRESOL)

! begin_doc_block
! ## `GPNORM_TRANSTL`
!
! ### Signature
!
! ```f90
! SUBROUTINE GPNORM_TRANSTL(PGP, KFIELDS, KPROMA, PAVE, KRESOL)
! ```
!
! ### Purpose
!
! This subroutine is the "tangent linear" version of `GPNORM_TRANS`. It is actually identical, but
! has a simpler interface as the `PMIN` and `PMAX` arguments are not used for tangent-linear
! applications.
!
! ### `INTENT(IN)` arguments
!
! - `REAL(KIND=JPRB), INTENT(IN) :: PGP(:,:,:)`  
!   Input grid point array, distributed across MPI tasks as usual.  
!   Dimensions: (NPROMA, number of fields, number of NPROMA blocks)
! - `INTEGER(KIND=JPIM), INTENT(IN) :: KFIELDS`  
!   Number of input fields on which to compute statistics.
! - `INTEGER(KIND=JPIM), INTENT(IN) :: KPROMA`  
!   Blocking factor for grid point input.
!
! ### `OPTIONAL, INTENT(IN)` arguments
!
! - `INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KRESOL`  
!   Resolution handle returned by original call to `SETUP_TRANS`.  
!   *Default*: `1` (i.e. first resolution handle)
!
! ### `INTENT(OUT)` arguments
!
! - `REAL(KIND=JPRB), INTENT(OUT) :: PAVE(:)`  
!   The global average of each field. Only task 1 will have data.  
!   Dimensions: (`KFIELDS`)
!
! end_doc_block

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

! Declaration of arguments
  
REAL(KIND=JPRB),INTENT(IN)    :: PGP(:,:,:)
REAL(KIND=JPRB),INTENT(OUT)   :: PAVE(:)
INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS
INTEGER(KIND=JPIM),INTENT(IN) :: KPROMA
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KRESOL

END SUBROUTINE GPNORM_TRANSTL
END INTERFACE
