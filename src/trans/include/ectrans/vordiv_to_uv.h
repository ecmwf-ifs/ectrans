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
SUBROUTINE VORDIV_TO_UV(PSPVOR,PSPDIV,PSPU,PSPV,KSMAX,KVSETUV)

! begin_doc_block
! ## `VORDIV_TO_UV`
!
! ### Signature
!
! ```f90
! SUBROUTINE VORDIV_TO_UV(PSPVOR, PSPDIV, PSPU, PSPV, KSMAX, KVSETUV)
! ```
!
! ### Purpose
!
! This subroutine converts spectral space vorticity and divergence to spectral space u and v wind
! components.
!
! Note that this is a "special" subroutine: it can be called in isolation without initialisation by
! `SETUP_TRANS0` and `SETUP_TRANS`. This is why it is necessary to pass in `KSMAX` explicitly.
!
! ### `INTENT(IN)` arguments
!
! - `REAL(KIND=JPRB), INTENT(IN) :: PSPVOR(:,:)`  
!   Spectral space vorticity.  
!   Dimensions: (vertical levels, spectral coefficients).
! - `REAL(KIND=JPRB), INTENT(IN) :: PSPDIV(:,:)`  
!   Spectral space divergence.  
!   Dimensions: (vertical levels, spectral coefficients).
! - `INTEGER(KIND=JPIM), INTENT(IN) :: KSMAX`  
!   Spectral truncation.
!
! ### `INTENT(OUT)` arguments
!
! - `REAL(KIND=JPRB), INTENT(OUT) :: PSPU(:,:)`  
!   Spectral space U (zonal) wind.  
!   Dimensions: (vertical levels, spectral coefficients).
! - `REAL(KIND=JPRB), INTENT(OUT) :: PSPV(:,:)`  
!   Spectral space V (meridional) wind.  
!   Dimensions: (vertical levels, spectral coefficients).
!
! ### `OPTIONAL, INTENT(IN)` arguments
!
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KVSETUV(:)`  
!   Array which maps each vertical level of the input vorticity/divergence to its corresponding  
!   member of the V-set.
! end_doc_block

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

! Declaration of arguments

REAL(KIND=JPRB), INTENT(IN) :: PSPVOR(:,:)
REAL(KIND=JPRB), INTENT(IN) :: PSPDIV(:,:)
REAL(KIND=JPRB), INTENT(OUT) :: PSPU(:,:)
REAL(KIND=JPRB), INTENT(OUT) :: PSPV(:,:)
INTEGER(KIND=JPIM) , INTENT(IN) :: KSMAX
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETUV(:)

END SUBROUTINE VORDIV_TO_UV
END INTERFACE
