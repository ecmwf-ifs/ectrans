! (C) Copyright 2026- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE BACKENDS_MOD

USE EC_PARKIND, ONLY: JPIM

IMPLICIT NONE

! Integer parameters to identify the backend type for which a deallocator is defined
INTEGER(KIND=JPIM), PARAMETER :: JP_BACKEND_CPU_SP = 1
INTEGER(KIND=JPIM), PARAMETER :: JP_BACKEND_CPU_DP = 2
INTEGER(KIND=JPIM), PARAMETER :: JP_BACKEND_GPU_SP = 3
INTEGER(KIND=JPIM), PARAMETER :: JP_BACKEND_GPU_DP = 4

! Currently there are 4 backends supported for global transforms
! Note: etrans currently still uses its own ETRANS_END finalise subroutine
INTEGER(KIND=JPIM), PARAMETER :: JP_N_BACKENDS = 4

! All final teardown subroutines must look like this
ABSTRACT INTERFACE
  SUBROUTINE BACKEND_TEARDOWN
  END SUBROUTINE BACKEND_TEARDOWN
END INTERFACE

TYPE :: BACKEND
  LOGICAL :: ENABLED = .FALSE.
  PROCEDURE(BACKEND_TEARDOWN), POINTER, NOPASS :: TEARDOWN => NULL()
END TYPE

TYPE(BACKEND) :: Y_BACKENDS(JP_N_BACKENDS)

END MODULE BACKENDS_MOD
