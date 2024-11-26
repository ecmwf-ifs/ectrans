! (C) Copyright 2022- NVIDIA.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE OPENACC_EXT_TYPE_MOD
  USE ISO_C_BINDING, ONLY: C_SIZE_T

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: EXT_ACC_ARR_DESC

  TYPE EXT_ACC_ARR_DESC
    INTEGER(C_SIZE_T) :: PTR, SZ
  END TYPE
END MODULE OPENACC_EXT_TYPE_MOD
