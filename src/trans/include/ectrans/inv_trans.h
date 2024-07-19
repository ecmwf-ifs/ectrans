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
SUBROUTINE INV_TRANS(PSPVOR,PSPDIV,PSPSCALAR,PSPSC3A,PSPSC3B,PSPSC2,&
 & FSPGL_PROC,&
 & LDSCDERS,LDVORGP,LDDIVGP,LDUVDER,LDLATLON,KPROMA,KVSETUV,KVSETSC,KRESOL,&
 & KVSETSC3A,KVSETSC3B,KVSETSC2,&
 & PGP,PGPUV,PGP3A,PGP3B,PGP2)

! begin_doc_block
! ## `INV_TRANS`
!
! ### Signature
!
! ```f90
! SUBROUTINE INV_TRANS(PSPVOR, PSPDIV, PSPSCALAR, PSPSC3A, PSPSC3B, PSPSC2, &
!   &                  FSPGL_PROC, LDSCDERS, LDVORGP, LDDIVGP, LDUVDER, LDLATLON, KPROMA, KVSETUV, &
!   &                  KVSETSC, KRESOL, KVSETSC3A, KVSETSC3B, KVSETSC2, PGP, PGPUV, PGP3A, PGP3B, &
!   &                  PGP2)
! ```
!
! ### Purpose
!
! This subroutine performs an inverse spectral transform.
! `INV_TRANS` supports two modes of passing arrays (which we call "call modes" below):
!
! 1. Pass in `PSPVOR`, `PSPDIV`, `PSPSCALAR`, receive `PGP`.
! 2. Pass in `PSPVOR`, `PSPDIV`, `PSPSC3A`, `PSPSC3B`, `PSPSC2`, receive `PGPUV`, `PGP3A`,  `PGP3B`, `PGP2`
!
! Although call mode 1 is simpler, it can often entail having to create array copies in order to
! subsequently operate on the different elements of `PGP`. Hence, with call mode 2, the output
! arrays are already categorized according to whether they are wind related, 3D scalar, or 2D
! scalar.
!
! Note that `PSPSC3A`/`PGP3A` and `PSPSC3B`/`PGP3B` are essentially identical in type and size.
! There are cases where it is useful to have the ability to pass two sets of scalar fields to be
! transformed independently.
!
! ### `OPTIONAL, INTENT(IN)` arguments
!
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PSPVOR(:,:)`  
!   Spectral space vorticity.  
!   Dimensions: (vertical levels, spectral coefficients).
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PSPDIV(:,:)`  
!   Spectral space divergence.  
!   Dimensions: (vertical levels, spectral coefficients).
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PSPSCALAR(:,:)`  
!   Spectral space scalar fields (call mode 1).  
!   Dimensions: (total number of 2D scalar fields, spectral coefficients).
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PSPSC3A(:,:,:)`  
!   Spectral space 3D scalar fields (call mode 2) corresponding to `PGP3A`.  
!   Dimensions: (vertical levels, spectral coefficients, number of 3D scalar fields).
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PSPSC3B(:,:,:)`  
!   Spectral space 3D scalar fields (call mode 2) corresponding to `PGP3B`.  
!   Dimensions: (vertical levels, spectral coefficients, number of 3D scalar fields).
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PSPSC2(:,:)`  
!   Spectral space 2D scalar fields (call mode 2) corresponding to `PGP2`.  
!   Dimensions: (number of 2D scalar fields, spectral coefficients).
! - `LOGICAL, OPTIONAL, INTENT(IN) :: LDSCDERS`  
!   Indicates if derivatives of scalar variables are required.  
!   *Default*: `.FALSE.`
! - `LOGICAL, OPTIONAL, INTENT(IN) :: LDVORGP`  
!   Indicates if grid point space vorticity is required.  
!   *Default*: `.FALSE.`
! - `LOGICAL, OPTIONAL, INTENT(IN) :: LDDIVGP`  
!   Indicates if grid point space divergence is required.  
!   *Default*: `.FALSE.`
! - `LOGICAL, OPTIONAL, INTENT(IN) :: LDUVDER`  
!   Indicates of East-West derivatives of U and V are required.  
!   *Default*: `.FALSE.`
! - `LOGICAL, OPTIONAL, INTENT(IN) :: LDLATLON`  
!   Indicates if regular latitude-longitude output is required.  
!   *Default*: `.FALSE.`
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KPROMA`  
!   Required blocking factor for grid point output.  
!   *Default*: `D%NGPTOT`
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KVSETUV(:)`  
!   Array which maps each vertical level of the input vorticity/divergence to its corresponding  
!   member of the V-set.
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KVSETSC(:)`  
!   Array which maps each of the input scalar fields to its corresponding member  of the V-set  
!   (call mode 1).
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KVSETSC3A(:)`  
!   Array which maps each vertical level of the input 3D scalar fields to its corresponding  
!   member of the V-set (call mode 2).
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KVSETSC3B(:)`  
!   Array which maps each vertical level of the input 3D scalar fields to its corresponding  
!   member of the V-set (call mode 2).
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KVSETSC2(:)`  
!   Array which maps each 2D scalar field to its corresponding member of the V-set (call mode 2).
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KRESOL`  
!   Resolution handle returned by original call to `SETUP_TRANS`.
!
! ### `OPTIONAL, EXTERNAL` arguments
! - `OPTIONAL, EXTERNAL :: FSPGL_PROC`  
!   External procedure to be executed in Fourier space before transposition.
!
! ### `OPTIONAL, INTENT(OUT)` arguments
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PGP(:,:,:)`  
!   Array containing all grid point fields (call mode 1).  
!   Dimensions: (NPROMA, number of grid point fields, number of NPROMA blocks).  
!   The ordering of grid point fields in this array is as follows:  
!     1. vorticity (if `PSPVOR` is present and `LDVORGP`)
!     2. divergence (if `PSPDIV` is present and `LDDIVGP`)
!     3. U wind (if `PSPVOR` and `PSPDIV` are present)
!     4. V wind (if `PSPVOR` and `PSPDIV` are present)
!     5. scalar fields (if `PSPSCALAR` is present)
!     6. North-South derivative of scalar fields (if `PSPSCALAR` is present and `LDSCDERS`)
!     7. East-West derivative of U wind (if `PSPVOR` and `PSPDIV` are present and `LDUVDER`)
!     8. East-West derivative of V wind (if `PSPVOR` and `PSPDIV` are present and `LDUVDER`)
!     9. East-West derivative of scalar fields (if `PSPSCALAR` is present and `LDSCDERS`)
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PGPUV(:,:,:,:)`  
!   Array containing grid point fields relating to wind (call mode 2).  
!   Dimensions: (NPROMA, vertical levels, number of wind-related fields, number of NPROMA blocks).  
!   The ordering of grid point fields in this array is as follows:  
!     1. vorticity (if `PSPVOR` is present and `LDVORGP`)
!     2. divergence (if `PSPDIV` is present and `LDDIVGP`)
!     3. U wind (if `PSPVOR` and `PSPDIV` are present)
!     4. V wind (if `PSPVOR` and `PSPDIV` are present)
!     5. East-West derivative of U wind (if `PSPVOR` and `PSPDIV` are present and `LDUVDER`)
!     6. East-West derivative of V wind (if `PSPVOR` and `PSPDIV` are present and `LDUVDER`)
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PGP3A(:,:,:,:)`  
!   Array containing 3D scalar grid point fields corresponding to `PSPSC3A`.  
!   Dimensions: (NPROMA, vertical levels, number of 3D scalar fields, number of NPROMA blocks).  
!   The ordering of grid point fields in this array is as follows:  
!     1. scalar fields
!     2. North-South derivatives of scalar fields (if `LDSCDERS`)
!     3. East-West derivatives of scalar fields (if `LDSCDERS`)
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PGP3B(:,:,:,:)`  
!   Array containing 3D scalar grid point fields corresponding to `PSPSC3B`.  
!   Dimensions: (NPROMA, vertical levels, number of 3D scalar fields, number of NPROMA blocks).  
!   The ordering of grid point fields in this array is as follows:  
!     1. scalar fields
!     2. North-South derivatives of scalar fields (if `LDSCDERS`)
!     3. East-West derivatives of scalar fields (if `LDSCDERS`)
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(OUT) :: PGP2(:,:,:)`  
!   Array containing 2D scalar grid point fields corresponding to `PSPSC2`.  
!   Dimensions: (NPROMA, number of 2D scalar fields, number of NPROMA blocks).  
!   The ordering of grid point fields in this array is as follows:  
!     1. scalar fields
!     2. North-South derivatives of scalar fields (if `LDSCDERS`)
!     3. East-West derivatives of scalar fields (if `LDSCDERS`)
! end_doc_block

USE PARKIND1  ,ONLY : JPIM     ,JPRB


IMPLICIT NONE

! Declaration of arguments

REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN) :: PSPVOR(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN) :: PSPDIV(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN) :: PSPSCALAR(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN) :: PSPSC3A(:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN) :: PSPSC3B(:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN) :: PSPSC2(:,:)
LOGICAL   ,OPTIONAL, INTENT(IN) :: LDSCDERS
LOGICAL   ,OPTIONAL, INTENT(IN) :: LDVORGP
LOGICAL   ,OPTIONAL, INTENT(IN) :: LDDIVGP
LOGICAL   ,OPTIONAL, INTENT(IN) :: LDUVDER
LOGICAL   ,OPTIONAL, INTENT(IN) :: LDLATLON
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KPROMA
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETUV(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC3A(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC3B(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC2(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KRESOL
EXTERNAL  FSPGL_PROC
OPTIONAL  FSPGL_PROC
REAL(KIND=JPRB),OPTIONAL    ,INTENT(OUT) :: PGP(:,:,:)
REAL(KIND=JPRB),OPTIONAL    ,INTENT(OUT) :: PGPUV(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    ,INTENT(OUT) :: PGP3A(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    ,INTENT(OUT) :: PGP3B(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    ,INTENT(OUT) :: PGP2(:,:,:)


END SUBROUTINE INV_TRANS

END INTERFACE
