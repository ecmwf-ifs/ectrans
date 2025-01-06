PROGRAM TEST_PROGRAM

USE PARKIND1, ONLY: JPIM, JPRB
USE MPL_MODULE

IMPLICIT NONE

! Spectral truncation
INTEGER(JPIM) :: TRUNC = 79
INTEGER(JPIM) :: verbosity = 1

! Arrays for storing our field in spectral space and grid point space
REAL(KIND=JPRB), ALLOCATABLE :: SPECTRAL_FIELD(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: GRID_POINT_FIELD(:,:,:)

! Dimensions of our arrays in spectral space and grid point space
INTEGER(KIND=JPIM) :: NSPEC2
INTEGER(KIND=JPIM) :: NGPTOT

#include "setup_trans0.h"
#include "setup_trans.h"
#include "trans_inq.h"
#include "inv_trans.h"

CALL MPL_INIT(ldinfo=(verbosity>=1))

CALL DR_HOOK_INIT()

! Initialise ecTrans (resolution-agnostic aspects)
CALL SETUP_TRANS0(LDMPOFF=.TRUE., KPRINTLEV=2)

! Initialise ecTrans (resolution-specific aspects)
CALL SETUP_TRANS(KSMAX=TRUNC, KDGL=2 * (TRUNC + 1))

! Inquire about the dimensions in spectral space and grid point space
CALL TRANS_INQ(KSPEC2=NSPEC2, KGPTOT=NGPTOT)

! Allocate our work arrays
ALLOCATE(SPECTRAL_FIELD(1,NSPEC2))
ALLOCATE(GRID_POINT_FIELD(NGPTOT,1,1))

! Initialise our spectral field array
SPECTRAL_FIELD(:,:) = 0.0_JPRB

! Perform an inverse transform
CALL INV_TRANS(PSPSCALAR=SPECTRAL_FIELD, PGP=GRID_POINT_FIELD)

CALL MPL_END(ldmeminfo=.false.)

END PROGRAM TEST_PROGRAM
