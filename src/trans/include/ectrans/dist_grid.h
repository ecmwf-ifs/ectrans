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
SUBROUTINE DIST_GRID(PGPG,KPROMA,KFDISTG,KFROM,KRESOL,PGP,KSORT)

! begin_doc_block
! ## `DIST_GRID`
!
! ### Signature
!
! ```f90
! SUBROUTINE DIST_GRID(PGPG, KPROMA, KFDISTG, KFROM, KRESOL, PGP, KSORT)
! ```
!
! ### Purpose
! The subroutine distributes one or more fields resident entirely on a single MPI task among all
! other tasks according to the grid point decomposition used in ecTrans.
!
! ### `INTENT(IN)` arguments
!
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KFDISTG`  
!   The number of fields to be distributed.
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KFROM`  
!   Array which, for each field to be distributed, which MPI task will be sending the field.
!
! ### `OPTIONAL, INTENT(IN)` arguments
!
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PGPG`  
!   Array containing fields to be distributed.  
!   Dimensions: (number of global grid points, number of fields on this MPI task). 
!   Note that this is optional, because some tasks may only receive fields (and so they wouldn't
!   have any fields to offer through `PGPG`).
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KPROMA`  
!   "Blocking factor" used for grid point output.
!   *Default*: `D%NGPTOT`
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KRESOL`  
!   Resolution handle returned by original call to `SETUP_TRANS`.  
!   *Default*: `1` (i.e. first resolution handle)
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KSORT`  
!   Array allowing to rearrange fields in the output array. For each element, specify which  
!   element you want the field to end up in.
!
! ### `INTENT(OUT)` arguments
!
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PGP(:,:)`  
!   Array containing grid point output.  
!   Dimensions: (NPROMA, number of fields, number of NPROMA blocks).
! end_doc_block

USE PARKIND1  ,ONLY : JPIM     ,JPRB


IMPLICIT NONE

! Declaration of arguments

REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGPG(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KPROMA
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFDISTG
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFROM(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KRESOL
REAL(KIND=JPRB)             , INTENT(OUT) :: PGP(:,:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KSORT (:)


!     ------------------------------------------------------------------

END SUBROUTINE DIST_GRID

END INTERFACE
