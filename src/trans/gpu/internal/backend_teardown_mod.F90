! (C) Copyright 2026- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE BACKEND_TEARDOWN_MOD

USE GROWING_ALLOCATOR_MOD, ONLY: DESTROY_GROWING_ALLOCATOR
USE TPM_FLT,               ONLY: S, FLT_RESOL
USE TPM_FIELDS_GPU,        ONLY: FIELDS_GPU_RESOL, FG
USE TPM_TRANS,             ONLY: GROWING_ALLOCATION

CONTAINS
SUBROUTINE BACKEND_TEARDOWN

IMPLICIT NONE

CALL DESTROY_GROWING_ALLOCATOR(GROWING_ALLOCATION)

NULLIFY(S)
IF (ALLOCATED(FLT_RESOL)) DEALLOCATE(FLT_RESOL)

NULLIFY(FG)
IF( ALLOCATED(FIELDS_GPU_RESOL) ) DEALLOCATE(FIELDS_GPU_RESOL)

END SUBROUTINE BACKEND_TEARDOWN
END MODULE BACKEND_TEARDOWN_MOD
