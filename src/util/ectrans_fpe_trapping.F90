! (C) Copyright 2025- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE ECTRANS_FPE_TRAPPING

IMPLICIT NONE

PRIVATE
PUBLIC :: ECTRANS_ENABLE_FPE

CONTAINS

SUBROUTINE ECTRANS_ENABLE_FPE()
  INTERFACE
    SUBROUTINE C_ECTRANS_ENABLE_FPE() BIND(C,NAME="ectrans_enable_fpe")
    END SUBROUTINE C_ECTRANS_ENABLE_FPE
  END INTERFACE

  CALL C_ECTRANS_ENABLE_FPE()
END SUBROUTINE ECTRANS_ENABLE_FPE

END MODULE ECTRANS_FPE_TRAPPING
