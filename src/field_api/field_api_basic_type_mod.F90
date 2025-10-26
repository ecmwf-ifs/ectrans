! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE FIELD_API_BASIC_TYPE_MOD

USE FIELD_BASIC_MODULE, ONLY: FIELD_BASIC

IMPLICIT NONE

TYPE FIELD_BASIC_PTR
 ! A pointer to field API field with additional METADA
  CLASS (FIELD_BASIC), POINTER :: PTR      ! POINTER TO FIELD_BASIC from field API
  INTEGER, ALLOCATABLE         :: IVSET(:) ! b-set  for spectral fields
  CHARACTER(LEN=10)            :: NAME     ! Name
END TYPE

END MODULE FIELD_API_BASIC_TYPE_MOD

