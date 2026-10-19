! (C) Copyright 2026- ECMWF.
! (C) Copyright 2026- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE TPMALD_FIELDS

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

TYPE ALDFIELDS_TYPE

REAL(KIND=JPRB) ,POINTER :: RLEPINM(:) ! eigen-values of the inverse Laplace operator
END TYPE ALDFIELDS_TYPE

TYPE(ALDFIELDS_TYPE),ALLOCATABLE,TARGET :: ALDFIELDS_RESOL(:)
TYPE(ALDFIELDS_TYPE),POINTER     :: FALD

REAL(KIND=JPRB) ,ALLOCATABLE :: FALD_RLEPINM(:)

END MODULE TPMALD_FIELDS
