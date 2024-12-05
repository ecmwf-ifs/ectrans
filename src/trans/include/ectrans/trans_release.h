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
SUBROUTINE TRANS_RELEASE(KRESOL)

! begin_doc_block
! ## `TRANS_RELEASE`
!
! ### Signature
!
! ```f90
! SUBROUTINE TRANS_RELEASE(KRESOL)
! ```
!
! ### Purpose
!
! This subroutine releases (i.e. deallocates and nullifies) all arrays related to the given
! resolution tag `KRESOL`.
!
! ### `INTENT(IN)` arguments
!
! - `INTEGER(KIND=JPIM), INTENT(IN) :: KRESOL`  
!   The handle for the resolution you want to release.
! end_doc_block

USE PARKIND1 ,ONLY : JPIM
INTEGER(KIND=JPIM),INTENT(IN) :: KRESOL
END SUBROUTINE TRANS_RELEASE
END INTERFACE
