! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE TPMALD_GEO

! Module containing data describing plane projection grid.

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

TYPE ALDGEO_TYPE

! GEOGRAPHY 
  
REAL(KIND=JPRB) :: EYWN       ! Y-reso 
REAL(KIND=JPRB) :: EXWN       ! X-reso
END TYPE ALDGEO_TYPE

TYPE(ALDGEO_TYPE),ALLOCATABLE,TARGET :: ALDGEO_RESOL(:)
TYPE(ALDGEO_TYPE),POINTER     :: GALD

END MODULE TPMALD_GEO
