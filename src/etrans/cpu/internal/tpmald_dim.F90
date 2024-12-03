! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE TPMALD_DIM

! Module for dimensions.

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

TYPE ALDDIM_TYPE

! COLLOCATION GRID DIMENSIONS
  
INTEGER(KIND=JPIM) :: NDGLSUR       ! Number of rows of latitudes+...
INTEGER(KIND=JPIM) :: NMSMAX        ! Zonal truncation
INTEGER(KIND=JPIM) :: NDGUX         ! Number of rows in zone C+I
END TYPE ALDDIM_TYPE

TYPE(ALDDIM_TYPE),ALLOCATABLE,TARGET :: ALDDIM_RESOL(:)
TYPE(ALDDIM_TYPE),POINTER     :: RALD

END MODULE TPMALD_DIM
