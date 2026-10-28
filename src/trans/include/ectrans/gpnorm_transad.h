! (C) Copyright 2024- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

INTERFACE
!> This subroutine is the "adjoint" version of GPNORM_TRANSTL.
!!
!! @param KFIELDS Number of input fields on which to compute statistics.
!! @param KPROMA Blocking factor for grid point input.
!! @param KRESOL Resolution handle returned by original call to SETUP_TRANS.
!!   Default: 1 (i.e. first resolution handle)
!! @param PAVE Input global average of each field. Only task 1 should have data. On output this will
!!   be zeroed.
!!   Dimensions: (KFIELDS)
!! @param PGP Output grid point array, distributed across MPI tasks as usual.
!!   Dimensions: (NPROMA, number of fields, number of NPROMA blocks)
SUBROUTINE GPNORM_TRANSAD(PGP,KFIELDS,KPROMA,PAVE,KRESOL)

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

! Declaration of arguments

REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGP(:,:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAVE(:)
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIELDS
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
INTEGER(KIND=JPIM),OPTIONAL, INTENT(IN)  :: KRESOL

END SUBROUTINE GPNORM_TRANSAD
END INTERFACE
