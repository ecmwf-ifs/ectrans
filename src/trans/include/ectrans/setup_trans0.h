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
!> This subroutine initialises all resolution-agnostic aspects of the ecTrans library. It must be
!! called before any other ecTrans subroutines.
!!
!! @param KOUT Fortran unit number for listing output.
!!   Default: 6 (i.e. STDOUT)
!! @param KERR Unit number for error messages.
!!   Default: 0 (i.e. STDERR)
!! @param KPRINTLEV Verbosity of standard output (0: output disabled, 1: normal, 2: debug).
!!   Default: 0
!! @param KMAX_RESOL Maximum number of different resolutions for which spectral transforms will be
!!   performed. This is required in advance so ecTrans knows how to correctly size the work arrays.
!!   Default: 1
!! @param KPROMATR Batch size for splitting of vertical field set. This allows one to transform one
!!   batch of vertical levels at a time rather than transforming all in one go. A value of 0
!!   disables this feature.
!!   Default: 0
!! @param KPRGPNS Splitting level in North-South direction in grid-point space.
!!   Default: 1
!! @param KPRGPEW Splitting level in East-West direction in grid-point space.
!!   Default: 1
!! @param KPRTRW Splitting level in wave direction in spectral space.
!!   Default: 1
!! @param LDMPOFF Switch off message passing.
!!   Default: .FALSE.
!! @param LDSYNC_TRANS Switch to activate barriers in M->L and L->M transposition routines.
!!   Default: .FALSE.
!! @param KTRANS_SYNC_LEVEL Control parameter for synchronization and blocking in communication
!!   routines.
!!   Default: 0
!! @param LDEQ_REGIONS Enable EQ_REGIONS partitioning mode.
!!   Default: .FALSE.
!! @param LDALLOPERM Allocate certain arrays permanently.
!!   Default: .FALSE.
!! @param PRAD Radius of the planet (Earth) in metres.
!!   Default: 6371229.0
!! @param KOPT_MEMORY_TR Memory strategy for gridpoint transpositions (0: heap, 1: stack).
!!   Default: 0
!! @param K_REGIONS Number of regions returned by EQ_REGIONS algorithm.
!! @param K_REGIONS_NS Maximum number of North-South partitions, as determined by EQ_REGIONS
!!   algorithm.
!! @param K_REGIONS_EW Maximum number of East-West partitions, as determined by EQ_REGIONS
!!   algorithm.
SUBROUTINE SETUP_TRANS0(KOUT,KERR,KPRINTLEV,KMAX_RESOL,KPROMATR,&
&                       KPRGPNS,KPRGPEW,KPRTRW,&
&                       LDMPOFF,LDSYNC_TRANS,KTRANS_SYNC_LEVEL,&
&                       LDEQ_REGIONS,K_REGIONS_NS,K_REGIONS_EW,K_REGIONS,&
&                       PRAD,LDALLOPERM,KOPT_MEMORY_TR)

USE EC_PARKIND  ,ONLY : JPIM     ,JPRD

IMPLICIT NONE

INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN)  :: KOUT,KERR,KPRINTLEV,KMAX_RESOL,KPROMATR
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN)  :: KPRGPNS,KPRGPEW,KPRTRW
LOGICAL            ,OPTIONAL,INTENT(IN)  :: LDMPOFF
LOGICAL            ,OPTIONAL,INTENT(IN)  :: LDSYNC_TRANS
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN)  :: KTRANS_SYNC_LEVEL
LOGICAL            ,OPTIONAL,INTENT(IN)  :: LDEQ_REGIONS
LOGICAL            ,OPTIONAL,INTENT(IN)  :: LDALLOPERM
REAL(KIND=JPRD)    ,OPTIONAL,INTENT(IN)  :: PRAD
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN)  :: KOPT_MEMORY_TR
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(OUT) :: K_REGIONS(:)
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(OUT) :: K_REGIONS_NS
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(OUT) :: K_REGIONS_EW

END SUBROUTINE SETUP_TRANS0



END INTERFACE
