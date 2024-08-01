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
SUBROUTINE DIR_TRANS(PSPVOR,PSPDIV,PSPSCALAR,PSPSC3A,PSPSC3B,PSPSC2,&
& LDLATLON,KPROMA,KVSETUV,KVSETSC,KRESOL,KVSETSC3A,KVSETSC3B,KVSETSC2,&
& PGP,PGPUV,PGP3A,PGP3B,PGP2)

! begin_doc_block
! ## `DIR_TRANS`
!
! ### Signature
!
! ```f90
! SUBROUTINE DIR_TRANS(PSPVOR, PSPDIV, PSPSCALAR, PSPSC3A, PSPSC3B, PSPSC2, LDLATLON, KPROMA, &
!   &                  KVSETUV, KVSETSC, KRESOL, KVSETSC3A, KVSETSC3B, KVSETSC2, PGP, PGPUV, &
!   &                  PGP3A, PGP3B, PGP2)
! ```
!
! ### Purpose
!
! This subroutine performs a direct spectral transform.
! `DIR_TRANS` supports two modes of passing arrays (which we call "call modes" below):
!
! 1. Pass in `PGP`, receive `PSPVOR`, `PSPDIV`, `PSPSCALAR`.
! 2. Pass in `PGPUV`, `PGP3A`,  `PGP3B`, `PGP2`, receive `PSPVOR`, `PSPDIV`, `PSPSC3A`, `PSPSC3B`, `PSPSC2`.
!
! Although call mode 1 is simpler, it can often entail having to create array copies in order to
! prepare the different elements of `PGP`. Hence, with call mode 2, the input and output arrays are
! categorized according to whether they are wind related, 3D scalar, or 2D scalar.
!
! Note that `PSPSC3A`/`PGP3A` and `PSPSC3B`/`PGP3B` are essentially identical in type and size.
! There are cases where it is useful to have the ability to pass two sets of scalar fields to be
! transformed independently.
!
! ### `OPTIONAL, INTENT(IN)` arguments
!
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KPROMA`  
!   Required blocking factor for grid point output.  
!   *Default*: `D%NGPTOT`
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KVSETUV(:)`  
!   Array which maps each vertical level of the output vorticity/divergence to its corresponding  
!   member of the V-set.
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KVSETSC(:)`  
!   Array which maps each of the output scalar fields to its corresponding member of the V-set  
!   (call mode 1).
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KVSETSC3A(:)`  
!   Array which maps each vertical level of the output 3D scalar fields to its corresponding  
!   member of the V-set (call mode 2).
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KVSETSC3B(:)`  
!   Array which maps each vertical level of the output 3D scalar fields to its corresponding  
!   member of the V-set (call mode 2).
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KVSETSC2(:)`  
!   Array which maps each output 2D scalar field to its corresponding member of the V-set  
!   (call mode 2).
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KRESOL`  
!   Resolution handle returned by original call to `SETUP_TRANS`.
! - `LOGICAL, OPTIONAL, INTENT(IN) :: LDLATLON`  
!   **TODO: what is this?**
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PGP(:,:,:)`  
!   Array containing all grid point fields (call mode 1).  
!   Dimensions: (NPROMA, number of grid point fields, number of NPROMA blocks).  
!   The ordering of grid point fields in this array is as follows:
!     1. U wind
!     2. V wind
!     3. scalar fields 
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PGPUV(:,:,:,:)`  
!   Array containing wind-related grid point fields (call mode 2).  
!   Dimensions: (NPROMA, vertical levels, number of wind fields, number of NPROMA blocks).  
!   `DIR_TRANS` only operates on U/V winds, not vorticity/divergence. Hence, the third dimension  
!   of `PGPUV` must be equal to `2` (one for U and one for V).
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PGP3A(:,:,:,:)`  
!   Array containing 3D scalar grid point fields, corresponding to `PSPSC3A`.  
!   Dimensions: (NPROMA, vertical levels, number of scalar fields, number of NPROMA blocks).  
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PGP3B(:,:,:,:)`  
!   Array containing 3D scalar grid point fields, corresponding to `PSPSC3B`.  
!   Dimensions: (NPROMA, vertical levels, number of scalar fields, number of NPROMA blocks).  
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PGP2(:,:,:)`  
!   Array containing 2D scalar grid point fields, corresponding to `PSPSC2`.  
!   Dimensions: (NPROMA, number of scalar fields, number of NPROMA blocks).
!
! ### `OPTIONAL, INTENT(OUT)` arguments
!
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)`  
!   Spectral space vorticity.
!   Dimensions: (vertical levels, spectral coefficients).
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)`  
!   Spectral space divergence.
!   Dimensions: (vertical levels, spectral coefficients).
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)`  
!   Spectral space scalar fields (call mode 1).  
!   Dimensions: (total number of 2D scalar fields, spectral coefficients).
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PSPSC3A(:,:,:)`  
!   Spectral space 3D scalar fields (call mode 2) corresponding to `PGP3A`.  
!   Dimensions: (vertical levels, spectral coefficients, number of 3D scalar fields).
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PSPSC3B(:,:,:)`  
!   Spectral space 3D scalar fields (call mode 2) corresponding to `PGP3B`.  
!   Dimensions: (vertical levels, spectral coefficients, number of 3D scalar fields).
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PSPSC2(:,:)`  
!   Spectral space 2D scalar fields (call mode 2) corresponding to `PGP2`.
! end_doc_block

USE PARKIND1  ,ONLY : JPIM     ,JPRB


IMPLICIT NONE

! Declaration of arguments

REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSC3A(:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSC3B(:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSC2(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KPROMA
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETUV(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC3A(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC3B(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC2(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KRESOL
LOGICAL   ,OPTIONAL, INTENT(IN) :: LDLATLON

REAL(KIND=JPRB),OPTIONAL    ,INTENT(IN) :: PGP(:,:,:)
REAL(KIND=JPRB),OPTIONAL    ,INTENT(IN) :: PGPUV(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    ,INTENT(IN) :: PGP3A(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    ,INTENT(IN) :: PGP3B(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    ,INTENT(IN) :: PGP2(:,:,:)


END SUBROUTINE DIR_TRANS

END INTERFACE
