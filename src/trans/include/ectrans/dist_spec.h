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
SUBROUTINE DIST_SPEC(PSPECG,KFDISTG,KFROM,KVSET,KRESOL,PSPEC,&
 & LDIM1_IS_FLD,KSMAX,KSORT)

! begin_doc_block
! ## `DIST_SPEC`
!
! ### Signature
!
! ```f90
! SUBROUTINE DIST_SPEC(PSPECG, KFDISTG, KFROM, KVSET, KRESOL, PSPEC, LDIM1_IS_FLD, KSMAX, KSORT)
! ```
!
! ### Purpose
!
! This subroutine distributes a global spectral array among MPI tasks according to the specified
! distribution parameters.
!
! ### `INTENT(IN)` arguments
!
! - `REAL(KIND=JPRB), INTENT(IN) :: KFDISTG`  
!   The global number of fields to be distributed.
! - `INTEGER(KIND=JPIM), INTENT(IN) :: KFROM(:)`  
!   Array which, for each field to be distributed, which MPI task is responsible for distributing
!   the field.  
!   Dimensions: (KFDISTG).
!
! ### `OPTIONAL, INTENT(IN)` arguments
!
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PSPECG(:,:)`  
!   Global spectral array to be distributed.  
!   Dimensions: (number of fields, number of spectral coefficients) (flipped if `LDIM1_IS_FLD` is  
!  `.FALSE.`).  
!   Note that this is optional, because some tasks may only receive fields (and so they wouldn't  
!   have any fields to offer through `PSPECG`).
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KVSET(:)`  
!   Array which, for each field to be distributed, which "V-set" the field belongs to.  
!   Dimensions: (KFDISTG).
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KRESOL`  
!   Resolution handle returned by original call to `SETUP_TRANS`.  
!   *Default*: `1` (i.e. first resolution handle)
! - `LOGICAL, OPTIONAL, INTENT(IN) :: LDIM1_IS_FLD`  
!   If `.TRUE.`, the first dimension of `PSPECG` and `PSPEC` corresponds to fields and the second  
!   dimension corresponds to spectral coefficients. If `.FALSE.`, the first dimension  
!   corresponds to spectral coefficients and the second dimension corresponds to fields.  
!   *Default*: `.TRUE.`
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KSMAX`  
!   Maximum total wavenumber of the spectral transform. If not provided, it is inferred from the  
!   resolution handle `KRESOL`.  
!   *Default*: inferred from `KRESOL`
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KSORT(:)`  
!   Array which, for each field to be distributed, specifies the sorting order of the spectral  
!   coefficients.  
!   Dimensions: (KFDISTG).  
!   *Default*: no sorting
!
! ### `OPTIONAL, INTENT(OUT)` arguments
!
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PSPEC(:,:)`  
!   Local spectral array on each MPI task after distribution.  
!   Dimensions: (number of fields on this MPI task, number of spectral coefficients on this MPI  
!   task).
! end_doc_block

USE PARKIND1  ,ONLY : JPIM     ,JPRB


IMPLICIT NONE

! Declaration of arguments

REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PSPECG(:,:)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFDISTG
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFROM(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSET(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KRESOL
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPEC(:,:)
LOGICAL            ,OPTIONAL, INTENT(IN)  :: LDIM1_IS_FLD
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KSMAX
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KSORT (:)


!     ------------------------------------------------------------------

END SUBROUTINE DIST_SPEC

END INTERFACE
