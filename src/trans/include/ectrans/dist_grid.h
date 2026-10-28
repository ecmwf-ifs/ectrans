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
!> The subroutine distributes one or more global grid point fields resident across one or more MPI
!! tasks (each field is entirely resident on a single MPI task without decomposition) among all
!! other tasks, decomposing them along the way. It is the opposite of GATH_GRID.
!!
!! The figure below illustrates an example in which four fields are distributed equally across four
!! MPI tasks. Note that at the end, each task receives part of all four input fields. Note also that
!! the NPROMA-based blocking is introduced at the end of DIST_GRID.
!!
!! ![A schematic showing how DIST_GRID works](img/dist_grid.png)
!!
!! @param KFDISTG The number of fields to be distributed.
!! @param KFROM Array which, for each field to be distributed, which MPI task will be sending the
!!   field.
!!   Dimensions: (KFDISTG)
!! @param PGPG Array containing fields to be distributed.
!!   Dimensions: (number of global grid points, number of fields on this MPI task).
!!   Note that this is optional, because some tasks may only receive fields (and so they wouldn't
!!   have any fields to offer through PGPG).
!! @param KPROMA "Blocking factor" used for grid point output.
!!   Default: D%NGPTOT
!! @param KRESOL Resolution handle returned by original call to SETUP_TRANS.
!!   Default: 1 (i.e. first resolution handle)
!! @param KSORT Array allowing to rearrange fields in the output array. For each element, specify
!!   which element you want the field to end up in.
!!   Dimensions: (second dimension of PGP)
!!   Default: no sorting
!! @param PGP Array containing grid point output.
!!   Dimensions: (NPROMA, number of fields, number of NPROMA blocks).
SUBROUTINE DIST_GRID(PGPG,KPROMA,KFDISTG,KFROM,KRESOL,PGP,KSORT)

USE PARKIND1  ,ONLY : JPIM     ,JPRB


IMPLICIT NONE

! Declaration of arguments

REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGPG(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KPROMA
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFDISTG
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFROM(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KRESOL
REAL(KIND=JPRB)             , INTENT(OUT) :: PGP(:,:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KSORT (:)


!     ------------------------------------------------------------------

END SUBROUTINE DIST_GRID

END INTERFACE
