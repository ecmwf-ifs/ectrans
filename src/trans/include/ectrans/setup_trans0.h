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
SUBROUTINE SETUP_TRANS0(KOUT,KERR,KPRINTLEV,KMAX_RESOL,KPROMATR,&
&                       KPRGPNS,KPRGPEW,KPRTRW,KCOMBFLEN,&
&                       LDMPOFF,LDSYNC_TRANS,KTRANS_SYNC_LEVEL,&
&                       LDEQ_REGIONS,K_REGIONS_NS,K_REGIONS_EW,K_REGIONS,&
&                       PRAD,LDALLOPERM,KOPT_MEMORY_TR)

! begin_doc_block
! ## `SETUP_TRANS0`
!
! ### Signature
!
! ```f90
! SUBROUTINE SETUP_TRANS0(KOUT, KERR, KPRINTLEV, KMAX_RESOL, KPROMATR, KPRGPNS, KPRGPEW, &
!   &                     KPRTRW, KCOMBFLEN, LDMPOFF, LDSYNC_TRANS, KTRANS_SYNC_LEVEL, &
!   &                     LDEQ_REGIONS, K_REGIONS_NS, K_REGIONS_EW, K_REGIONS, PRAD, LDALLOPERM, &
!   &                     KOPT_MEMORY_TR)
! ```
!
! ### Purpose
!
! This subroutine initialises all resolution-agnostic aspects of the ecTrans library. It must be
! called before any other ecTrans subroutines.
!
! ### `OPTIONAL, INTENT(IN)` arguments
!
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KOUT`  
!   Fortran unit number for listing output.  
!   *Default*: `6` (i.e. STDOUT)
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KERR`  
!   Unit number for error messages.  
!   *Default*: `0` (i.e. STDERR)
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KPRINTLEV`  
!   Verbosity of standard output (0: output disabled, 1: normal, 2: debug).  
!   *Default*: `0`
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KMAX_RESOL`  
!   Maximum number of different resolutions for which spectral transforms will be performed. This is
!   required in advance so ecTrans knows how to correctly size the work arrays.  
!   *Default*: `1`
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KPROMATR`  
!   Batch size for splitting of vertical field set. This allows one to transform one batch of vertical  
!   levels at a time rather than transforming all in one go. A value of 0 disables this feature.  
!   *Default*: `0`  
!   **Potentially deprecatable.**
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KPRGPNS`  
!   Splitting level in North-South direction in grid-point space.  
!   *Default*: `1`
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KPRGPEW`  
!   Splitting level in East-West direction in grid-point space.  
!   *Default*: `1`
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KPRTRW`  
!   Splitting level in wave direction in spectral space.  
!   *Default*: `1`
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KCOMBFLEN`  
!   Size of communication buffer [1800000 (8bytes) ].  
!   **This argument is deprecated and will be removed in the next release.**
! - `LOGICAL, OPTIONAL, INTENT(IN) :: LDMPOFF`  
!   Switch off message passing.  
!   *Default*: `.FALSE.`  
!   **Potentially deprecatable.**
! - `LOGICAL, OPTIONAL, INTENT(IN) :: LDSYNC_TRANS`  
!   Switch to activate barriers in M->L and L->M transposition routines.  
!   *Default*: `.FALSE.`
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KTRANS_SYNC_LEVEL`  
!   Control parameter for synchronization and blocking in communication routines.  
!   *Default*: `0`
! - `LOGICAL, OPTIONAL, INTENT(IN) :: LDEQ_REGIONS`  
!   Enable EQ_REGIONS partitioning mode.  
!   *Default*: `.FALSE.`
! - `LOGICAL, OPTIONAL, INTENT(IN) :: LDALLOPERM`  
!   Allocate certain arrays permanently.  
!   *Default*: `.FALSE.`
! - `REAL(KIND=JPRD), OPTIONAL, INTENT(IN) :: PRAD`  
!   Radius of the planet (Earth) in metres.  
!   *Default*: `6371229.0`
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KOPT_MEMORY_TR`  
!   Memory strategy for gridpoint transpositions (0: heap, 1: stack).  
!   *Default*: `0`
!
! ### `OPTIONAL, INTENT(OUT)` arguments
!
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: K_REGIONS(:)`  
!   Number of regions returned by EQ_REGIONS algorithm.
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: K_REGIONS_NS`  
!   Maximum number of North-South partitions, as determined by EQ_REGIONS algorithm.
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: K_REGIONS_EW`  
!   Maximum number of East-West partitions, as determined by EQ_REGIONS algorithm.
! end_doc_block

USE EC_PARKIND  ,ONLY : JPIM     ,JPRD

IMPLICIT NONE

INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN)  :: KOUT,KERR,KPRINTLEV,KMAX_RESOL,KPROMATR
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN)  :: KPRGPNS,KPRGPEW,KPRTRW,KCOMBFLEN
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
