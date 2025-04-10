! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE TPMALD_DISTR

! Module for distributed memory environment.

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

TYPE ALDDISTR_TYPE

INTEGER(KIND=JPIM) ,POINTER :: NESM0(:)  ! Address in a spectral array of (m, n=m)
INTEGER(KIND=JPIM) ,POINTER :: NCPL2M(:) ! Number of complex Laplace coefficient for m given
INTEGER(KIND=JPIM) ,POINTER :: NPME(:)   ! Address for the Laplace operator and its inverse

END TYPE ALDDISTR_TYPE

TYPE(ALDDISTR_TYPE),ALLOCATABLE,TARGET :: ALDDISTR_RESOL(:)
TYPE(ALDDISTR_TYPE),POINTER     :: DALD

END MODULE TPMALD_DISTR

