SUBROUTINE SETUP_TRANS(KSMAX,KDGL,KLOEN,LDSPLIT,KAPSETS)

!**** *SETUP_TRANS* - Setup transform package for specific resolution

!     Purpose.
!     --------
!     To setup for making spectral transforms

!**   Interface.
!     ----------
!     CALL SETUP_TRANS(...)

!     Explicit arguments : KLOEN,LDSPLIT,KAPSETS are optional arguments
!     -------------------- 
!     KSMAX - spectral truncation required
!     KDGL  - number of gaussian latitudes
!     KLOEN(:) - number of points on each Gaussian latitude
!     LDSPLIT - true if split latitudes in grid-point space
!     KAPSETS - Number of apple sets in the distribution
!
!     KSMAX,KDGL and KLOEN are GLOBAL varaibles desribing the resolition
!     in spectral and grid-point space

!     LDSPLIT and KAPSETS describe the distribution among processors of
!     grid-point data and has no relevance if you are using a single processor
 
!     ------------------------------------------------------------------

#include "tsmbkind.h"

!ifndef INTERFACE

USE TPM_GEN
USE TPM_DIM
USE TPM_DISTR
USE TPM_GEOMETRY
USE TPM_FIELDS
USE TPM_FFT

USE SET_RESOL_MOD
USE SETUP_DIMS_MOD
USE SUMP_TRANS_MOD
USE SUMP_TRANS_PRELEG_MOD
USE SULEG_MOD
USE SETUP_GEOM_MOD
USE SUFFT_MOD

!endif INTERFACE

IMPLICIT NONE

! Dummy arguments

INTEGER_M ,INTENT(IN) :: KSMAX,KDGL
INTEGER_M ,OPTIONAL,INTENT(IN) :: KLOEN(:)
LOGICAL   ,OPTIONAL,INTENT(IN) :: LDSPLIT
INTEGER_M ,OPTIONAL,INTENT(IN) :: KAPSETS

!ifndef INTERFACE

! Local variables
INTEGER_M :: JGL,IERROR

LOGICAL :: LLP1,LLP2

!     ------------------------------------------------------------------


IF(MSETUP0 /= 1) THEN
  CALL ABOR1('SETUP_TRANS: SETUP_TRANS0 HAS TO BE CALLED BEFORE SETUP_TRANS')
ENDIF

LLP1 = NPRINTLEV>0
LLP2 = NPRINTLEV>1
IF(LLP1) WRITE(NOUT,*) '=== ENTER ROUTINE SETUP_TRANS ==='

! Allocate resolution dependent structures
IF(.NOT. ALLOCATED(DIM_RESOL)) THEN
  NDEF_RESOL = 1
  ALLOCATE(DIM_RESOL(NMAX_RESOL))
  ALLOCATE(FIELDS_RESOL(NMAX_RESOL))
  ALLOCATE(GEOM_RESOL(NMAX_RESOL))
  ALLOCATE(DISTR_RESOL(NMAX_RESOL))
  ALLOCATE(FFT_RESOL(NMAX_RESOL))
ELSE
  NDEF_RESOL = NDEF_RESOL+1
  IF(NDEF_RESOL > NMAX_RESOL) THEN
    CALL ABOR1('SETUP_TRANS:NDEF_RESOL > NMAX_RESOL')
  ENDIF
ENDIF

! Point at structures due to be initialized
CALL SET_RESOL(NDEF_RESOL)

IF(LLP1) WRITE(NOUT,*) '=== DEFINING RESOLUTION ',NCUR_RESOL



! Defaults for optional arguments


G%LREDUCED_GRID = .FALSE.
G%LINEAR_GRID = .FALSE.
D%LSPLIT = .FALSE.
D%NAPSETS = 0

! NON-OPTIONAL ARGUMENTS
R%NSMAX = KSMAX
R%NDGL  = KDGL
R%NDLON = 2*KDGL

! Optional arguments


ALLOCATE(G%NLOEN(R%NDGL))
IF(LLP2)WRITE(NOUT,9) 'NLOEN   ',SIZE(G%NLOEN   ),SHAPE(G%NLOEN   )
IF(PRESENT(KLOEN)) THEN
  DO JGL=1,R%NDGL
    IF(KLOEN(JGL) .NE. R%NDLON) THEN
      G%LREDUCED_GRID = .TRUE.
      EXIT
    ENDIF
  ENDDO
ENDIF
IF (G%LREDUCED_GRID) THEN
  G%NLOEN(:) = KLOEN(1:R%NDGL)
ELSE
  G%NLOEN(:) = R%NDLON
ENDIF

IF(PRESENT(LDSPLIT)) THEN
  D%LSPLIT = LDSPLIT
ENDIF
IF(PRESENT(KAPSETS)) THEN
  D%NAPSETS = KAPSETS
ENDIF

!Temporary?
R%NTMAX = R%NSMAX
IF(R%NSMAX > (R%NDLON+3)/3) THEN
  G%LINEAR_GRID = .TRUE.
ENDIF  

!     Setup resolution dependent structures
!     -------------------------------------

! Setup distribution independent dimensions
CALL SETUP_DIMS

! First part of setup of distributed environment
CALL SUMP_TRANS_PRELEG

! Compute Legandre polonomial and Gaussian Latitudes and Weights
CALL SULEG

! Compute arrays related to grid-point geometry
CALL SETUP_GEOM

! Second part of setup of distributed environment
CALL SUMP_TRANS

! Initialize Fast Fourier Transform package
CALL SUFFT

!     ------------------------------------------------------------------
9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

!endif INTERFACE

END SUBROUTINE SETUP_TRANS


