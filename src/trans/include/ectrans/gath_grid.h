! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

INTERFACE
!> The subroutine gathers one or more fields distributed across one or more MPI tasks and sends them
!! to a particular MPI task. It is the opposite of DIST_GRID.
!!
!! The figure below illustrates an example in which four fields which are distributed equally across
!! four MPI tasks are gathered so that each task possesses one field in its entirety. The NPROMA
!! blocking is removed along the way so that at the end, each field only has a single horizontal
!! dimension
!!
!! ![A schematic showing how DIST_GRID works](img/gath_grid.png)
!!
!! @param KFGATHG The number of fields to be gathered.
!! @param KTO Array which, for each field to be gathered, which MPI task will be receiving the
!!   field.
!!   Dimensions: (KFGATHG)
!! @param PGP Array containing grid point output.
!!   Dimensions: (NPROMA, number of fields, number of NPROMA blocks).
!! @param KPROMA "Blocking factor" used for grid point output.
!!   Default: D%NGPTOT
!! @param KRESOL Resolution handle returned by original call to SETUP_TRANS.
!!   Default: 1 (i.e. first resolution handle)
!! @param PGPG Array containing gathered fields.
!!   Dimensions: (number of global grid points, number of fields on this MPI task).
!!   Note that this is optional, because not all tasks may be receiving fields.
SUBROUTINE GATH_GRID(PGPG,KPROMA,KFGATHG,KTO,KRESOL,PGP)

USE PARKIND1  ,ONLY : JPIM     ,JPRB


IMPLICIT NONE

! Declaration of arguments

REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PGPG(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KPROMA
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFGATHG
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KTO(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KRESOL
REAL(KIND=JPRB)             , INTENT(IN)  :: PGP(:,:,:)


!     ------------------------------------------------------------------

END SUBROUTINE GATH_GRID

END INTERFACE
