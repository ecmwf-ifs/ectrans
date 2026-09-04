! (C) Copyright 2026- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE RESOLS_MOD

USE EC_PARKIND, ONLY: JPIM

IMPLICIT NONE

! All deallocators must look like this
ABSTRACT INTERFACE
  SUBROUTINE DEALLOC_RESOL(KRESOL)
    USE EC_PARKIND, ONLY: JPIM
    INTEGER(KIND=JPIM), INTENT(IN) :: KRESOL
  END SUBROUTINE DEALLOC_RESOL
END INTERFACE

TYPE :: RESOL
  INTEGER(KIND=JPIM) :: BACKEND ! Backend ID (see BACKENDS_MOD)
  PROCEDURE(DEALLOC_RESOL), POINTER, NOPASS :: DESTROY => NULL()
END TYPE

TYPE(RESOL), ALLOCATABLE :: Y_RESOLS(:)

END MODULE RESOLS_MOD
