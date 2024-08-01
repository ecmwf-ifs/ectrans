! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

INTERFACE
SUBROUTINE GET_CURRENT(KRESOL,LDLAM)

! begin_doc_block
! ## `GET_CURRENT`
!
! ### Signature
!
! ```f90
! SUBROUTINE GET_CURRENT(KRESOL,LDLAM)
! ```
!
! ### Purpose
!
! This subroutine allows one to retrieve the current resolution handle, `NCUR_RESOL`, and the LAM
! parameter, `LDLAM`, which indicates if the local area version of ecTrans is being used or not.
!
! @note
! `LDLAM` will always be `.FALSE.` for the moment. ecTrans currently only supports global spectral
! transforms, though local (bifourier) transforms may be supported in future.
! @endnote
!
! ### `OPTIONAL, INTENT(OUT)` arguments
!
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KRESOL`  
!   The current default resolution handle.
! - `LOGICAL, OPTIONAL, INTENT(OUT) :: LDLAM`  
!   Whether the local area version of ecTrans is being used or not.
! end_doc_block

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

INTEGER(KIND=JPIM)  ,OPTIONAL,INTENT(OUT)  :: KRESOL
LOGICAL             ,OPTIONAL,INTENT(OUT)  :: LDLAM

END SUBROUTINE GET_CURRENT
END INTERFACE
