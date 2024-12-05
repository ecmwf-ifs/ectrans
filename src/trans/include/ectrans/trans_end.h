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
SUBROUTINE TRANS_END(CDMODE)

! begin_doc_block
! ## `TRANS_END`
!
! ### Signature
!
! ```f90
! SUBROUTINE TRANS_END(CDMODE)
! ```
!
! ### Purpose
!
! This subroutine terminates the transform package by releasing (i.e. deallocating and nullifying)
! all allocated arrays.
!
! ### `OPTIONAL, INTENT(IN)` arguments
!
! - `CHARACTER*5, OPTIONAL, INTENT(IN) :: CDMODE`  
!   A string parameter indicating which mode to finalise with. Either `'FINAL'` to finalise  
!   everything or `'INTER'` just to deallocate `N_REGIONS` and `NPRCIDS`.  
!   *Default*: `'FINAL'`
! end_doc_block

USE PARKIND1  ,ONLY : JPIM     ,JPRB
IMPLICIT NONE
CHARACTER*5, OPTIONAL, INTENT(IN) :: CDMODE


END SUBROUTINE TRANS_END
END INTERFACE
