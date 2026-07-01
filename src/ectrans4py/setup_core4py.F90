MODULE SETUP_CORE4PY_MOD
! Shared core for the ectrans4py setup, so the resolution-independent SETUP_TRANS0
! is called from a single place. Both the serial init (SPEC_SETUP4PY, LDMPOFF) and
! the distributed init (SETUP_TRANS0_4PY, processor grid) go through SETUP_TRANS0_CORE.
!
! The parallelism of ecTrans is fixed at SETUP_TRANS0 time (LDMPOFF and the
! processor grid KPRGPNS x KPRGPEW x KPRTRW), which is why the two setups keep
! separate Python entry points; this only removes the duplicated SETUP_TRANS0
! marshalling between them.
IMPLICIT NONE
PRIVATE
PUBLIC :: SETUP_TRANS0_CORE

CONTAINS

SUBROUTINE SETUP_TRANS0_CORE(KMAX_RESOL, LDMPOFF, KPRGPNS, KPRGPEW, KPRTRW, &
                           & LDEQ_REGIONS, KREGIONS_NS, KREGIONS_EW)
! Resolution-independent transform setup (IFS SUTRANS0 'SETUP' branch).
!    KMAX_RESOL : maximum number of resolutions
!    LDMPOFF    : .TRUE. => message passing off (serial, single task)
!    KPRGPNS/KPRGPEW/KPRTRW : grid-point (N-S x E-W) and wave-set processor grid
!                             (pass 1,1,1 for the serial case)
!    LDEQ_REGIONS           : equal-regions partitioning
!    KREGIONS_NS/KREGIONS_EW : returned region counts
USE PARKIND1, ONLY: JPIM
IMPLICIT NONE
INTEGER(KIND=JPIM), INTENT(IN)  :: KMAX_RESOL, KPRGPNS, KPRGPEW, KPRTRW
LOGICAL,            INTENT(IN)  :: LDMPOFF, LDEQ_REGIONS
INTEGER(KIND=JPIM), INTENT(OUT) :: KREGIONS_NS, KREGIONS_EW
INTEGER(KIND=JPIM) :: IREGIONS(KPRGPNS*KPRGPEW+2)
#include "setup_trans0.h"
IREGIONS(:) = 0
KREGIONS_NS = 0
KREGIONS_EW = 0
CALL SETUP_TRANS0(KPRGPNS=KPRGPNS, KPRGPEW=KPRGPEW, KPRTRW=KPRTRW, &
 & KMAX_RESOL=KMAX_RESOL, LDEQ_REGIONS=LDEQ_REGIONS, LDMPOFF=LDMPOFF, KPRINTLEV=0, &
 & K_REGIONS=IREGIONS, K_REGIONS_NS=KREGIONS_NS, K_REGIONS_EW=KREGIONS_EW)
END SUBROUTINE SETUP_TRANS0_CORE

END MODULE SETUP_CORE4PY_MOD
