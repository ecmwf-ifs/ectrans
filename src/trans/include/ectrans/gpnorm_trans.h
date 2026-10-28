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
!> This subroutine calculates the global average (and optionally the global minimum and maximum) of
!! the input fields in grid point space. The results are communicated to task 1 only.
!!
!! @param PGP Input grid point array, distributed across MPI tasks as usual.
!!   Dimensions: (NPROMA, number of fields, number of NPROMA blocks)
!! @param KFIELDS Number of input fields on which to compute statistics.
!! @param KPROMA Blocking factor for grid point input.
!! @param LDAVE_ONLY Specify that only the average should be calculated before the aggregation to
!!   task 1. In this case PMIN and PMAX should already contain the minimum and maximum values of
!!   each field for each MPI task's local region, and afterwards they will contain the global
!!   maximum and minimum (again, only on task 1).
!! @param KRESOL Resolution handle returned by original call to SETUP_TRANS.
!!   Default: 1 (i.e. first resolution handle)
!! @param LPGP_ON_GPU When .TRUE., this argument indicates that input grid point data is already
!!   resident on device and no host->device update is needed.
!!   Default: .FALSE.
!! @param PAVE The global average of each field. Only task 1 will have data.
!!   Dimensions: (KFIELDS)
!! @param PMIN The global minimum of each field. Only task 1 will have data.
!!   Dimensions: (KFIELDS)
!! @param PMAX The global maximum of each field. Only task 1 will have data.
!!   Dimensions: (KFIELDS)
SUBROUTINE GPNORM_TRANS(PGP,KFIELDS,KPROMA,PAVE,PMIN,PMAX,LDAVE_ONLY,KRESOL,LPGP_ON_GPU)

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

! Declaration of arguments
  
REAL(KIND=JPRB),INTENT(IN)    :: PGP(:,:,:)
REAL(KIND=JPRB),INTENT(OUT)   :: PAVE(:)
REAL(KIND=JPRB),INTENT(INOUT) :: PMIN(:)
REAL(KIND=JPRB),INTENT(INOUT) :: PMAX(:)
INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS
INTEGER(KIND=JPIM),INTENT(IN) :: KPROMA
LOGICAL,INTENT(IN)            :: LDAVE_ONLY
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KRESOL
LOGICAL            ,OPTIONAL, INTENT(IN)  :: LPGP_ON_GPU

END SUBROUTINE GPNORM_TRANS
END INTERFACE
