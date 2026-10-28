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
SUBROUTINE GATH_SPEC(PSPECG,KFGATHG,KTO,KVSET,KRESOL,PSPEC,LDIM1_IS_FLD,KSMAX,LDZA0IP)

! begin_doc_block
! ## `GATH_SPEC`
!
! ### Signature
!
! ```f90
! SUBROUTINE GATH_SPEC(PSPECG, KFGATHG, KTO, KVSET, KRESOL, PSPEC, LDIM1_IS_FLD, KSMAX, LDZA0IP)
! ```
!
! ### Purpose
!
! This subroutine gathers spectral fields decomposed among MPI tasks according to the specified
! distribution parameters. This subroutine is the opposite of `DIST_SPEC`.
!
! The figure below illustrates an example in which four fields distributed across four MPI
! tasks are collected so that each task has one global spectral field. In this case, `NPRTRV = 2`
! which means that each field is split among two MPI tasks to begin with, as shown.
!
! ![A schematic showing how GATH_SPEC works](img/gath_spec.png)
!
! ### `INTENT(IN)` arguments
!
! - `REAL(KIND=JPRB), INTENT(IN) :: KFGATHG`  
!   The global number of fields to be gathered.
! - `INTEGER(KIND=JPIM), INTENT(IN) :: KTO(:)`  
!   Array which, for each field to be gathered, which MPI task will receive the field.  
!   Dimensions: (KFGATHG).
!
! ### `OPTIONAL, INTENT(IN)` arguments
!
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PSPEC(:,:)`  
!   Array containing this MPI task's portion of the spectral fields to be gathered.  
!   Dimensions: (number of fields, number of spectral coefficients) (flipped if `LDIM1_IS_FLD` is  
!  `.FALSE.`).  
!   Note that this is optional, because some tasks may not actually own any spectral fields (and  
!   so they wouldn't have any fields to offer through `PSPEC`).
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KVSET(:)`  
!   Array which, for each field to be distributed, which "V-set" the field belongs to.  
!   Dimensions: (KFGATHG).
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
! - `LOGICAL, OPTIONAL, INTENT(IN) :: LDZA0IP`  
!   IF `.TRUE.`, the imaginary part of the coefficients corresponding to zonal wavenumber 0 are  
!   zeroed.  
!   *Default*: `.TRUE.`
!
! ### `OPTIONAL, INTENT(OUT)` arguments
!
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PSPECG(:,:)`  
!   Array containing gathered fields.  
!   Dimensions: (number of global spectral coefficients, number of fields on this MPI task).  
!   Note that this is optional, because not all tasks may be receiving fields.
! end_doc_block

USE PARKIND1  ,ONLY : JPIM     ,JPRB


IMPLICIT NONE

! Declaration of arguments

REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPECG(:,:)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFGATHG
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KTO(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSET(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KRESOL
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PSPEC(:,:)
LOGICAL            ,OPTIONAL, INTENT(IN)  :: LDIM1_IS_FLD
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KSMAX
LOGICAL            ,OPTIONAL, INTENT(IN)  :: LDZA0IP


!     ------------------------------------------------------------------

END SUBROUTINE GATH_SPEC

END INTERFACE
