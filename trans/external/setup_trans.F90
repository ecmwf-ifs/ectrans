SUBROUTINE SETUP_TRANS(KSMAX,KDGL,KLOEN,LDLINEAR_GRID,LDSPLIT,&
&KAPSETS,KTMAX,KRESOL,LDGRIDONLY)

!**** *SETUP_TRANS* - Setup transform package for specific resolution

!     Purpose.
!     --------
!     To setup for making spectral transforms. Each call to this routine
!     creates a new resolution up to a maximum of NMAX_RESOL set up in
!     SETUP_TRANS0. You need to call SETUP_TRANS0 before this routine can
!     be called.

!**   Interface.
!     ----------
!     CALL SETUP_TRANS(...)

!     Explicit arguments : KLOEN,LDLINEAR_GRID,LDSPLIT,KAPSETS are optional arguments
!     -------------------- 
!     KSMAX - spectral truncation required
!     KDGL  - number of Gaussian latitudes
!     KLOEN(:) - number of points on each Gaussian latitude [2*KDGL]
!     LDSPLIT - true if split latitudes in grid-point space [false]
!     LDLINEAR_GRID - true if linear grid
!     KAPSETS - Number of apple sets in the distribution [0]
!     KTMAX - truncation order for tendencies?
!     KRESOL - the resolution identifier
!     LDGRIDONLY - true if only grid space is required

!     KSMAX,KDGL,KTMAX and KLOEN are GLOBAL variables desribing the resolution
!     in spectral and grid-point space

!     LDSPLIT and KAPSETS describe the distribution among processors of
!     grid-point data and has no relevance if you are using a single processor
 
!     Method.
!     -------

!     Externals.  SET_RESOL   - set resolution
!     ----------  SETUP_DIMS  - setup distribution independent dimensions
!                 SUMP_TRANS_PRELEG - first part of setup of distr. environment
!                 SULEG - Compute Legandre polonomial and Gaussian 
!                         Latitudes and Weights
!                 SETUP_GEOM - Compute arrays related to grid-point geometry
!                 SUMP_TRANS - Second part of setup of distributed environment
!                 SUFFT - setup for FFT

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

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
USE ABORT_TRANS_MOD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

!endif INTERFACE

IMPLICIT NONE

! Dummy arguments

INTEGER(KIND=JPIM) ,INTENT(IN) :: KSMAX,KDGL
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN) :: KLOEN(:)
LOGICAL   ,OPTIONAL,INTENT(IN) :: LDLINEAR_GRID
LOGICAL   ,OPTIONAL,INTENT(IN) :: LDSPLIT
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN) :: KAPSETS
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN) :: KTMAX
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(OUT):: KRESOL
LOGICAL   ,OPTIONAL,INTENT(IN):: LDGRIDONLY

!ifndef INTERFACE

! Local variables
INTEGER(KIND=JPIM) :: JGL

LOGICAL :: LLP1,LLP2
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SETUP_TRANS',0,ZHOOK_HANDLE)

IF(MSETUP0 /= 1) THEN
  CALL ABORT_TRANS('SETUP_TRANS: SETUP_TRANS0 HAS TO BE CALLED BEFORE SETUP_TRANS')
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
    CALL ABORT_TRANS('SETUP_TRANS:NDEF_RESOL > NMAX_RESOL')
  ENDIF
ENDIF

IF (PRESENT(KRESOL)) THEN
  KRESOL=NDEF_RESOL
ENDIF

! Point at structures due to be initialized
CALL SET_RESOL(NDEF_RESOL)

IF(LLP1) WRITE(NOUT,*) '=== DEFINING RESOLUTION ',NCUR_RESOL



! Defaults for optional arguments


G%LREDUCED_GRID = .FALSE.
G%LINEAR_GRID = .FALSE.
D%LGRIDONLY = .FALSE.
D%LSPLIT = .FALSE.
D%NAPSETS = 0

! NON-OPTIONAL ARGUMENTS
R%NSMAX = KSMAX
R%NDGL  = KDGL
R%NDLON = 2*KDGL

IF (KDGL <= 0 .OR. MOD(KDGL,2) /= 0) THEN
  CALL ABORT_TRANS ('SETUP_TRANS: KDGL IS NOT A POSITIVE, EVEN NUMBER')
ENDIF

! Optional arguments

ALLOCATE(G%NLOEN(R%NDGL))
IF(LLP2)WRITE(NOUT,9) 'NLOEN   ',SIZE(G%NLOEN   ),SHAPE(G%NLOEN   )
IF(PRESENT(KLOEN)) THEN
  DO JGL=1,R%NDGL
    IF(KLOEN(JGL) /= R%NDLON) THEN
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

IF(PRESENT(KTMAX)) THEN
  R%NTMAX = KTMAX
ELSE
  R%NTMAX = R%NSMAX
ENDIF
IF(R%NTMAX /= R%NSMAX) THEN
  !This SHOULD work but I don't know how to test it /MH
  CALL ABORT_TRANS('SETUP_TRANS:R%NTMAX /= R%NSMAX HAS NOT BEEN VALIDATED')
ENDIF
!Temporary?
IF(PRESENT(LDLINEAR_GRID)) THEN
  G%LINEAR_GRID = LDLINEAR_GRID
ELSEIF(R%NSMAX > (R%NDLON+3)/3) THEN
  G%LINEAR_GRID = .TRUE.
ENDIF  

IF(PRESENT(LDGRIDONLY)) THEN
  D%LGRIDONLY=LDGRIDONLY
ENDIF

!     Setup resolution dependent structures
!     -------------------------------------

! Setup distribution independent dimensions
CALL SETUP_DIMS

! First part of setup of distributed environment
CALL SUMP_TRANS_PRELEG

! Compute Legandre polonomial and Gaussian Latitudes and Weights
CALL SULEG

CALL GSTATS(1802,0)

! Compute arrays related to grid-point geometry
CALL SETUP_GEOM
CALL GSTATS(1802,1)

! Second part of setup of distributed environment
CALL SUMP_TRANS
CALL GSTATS(1802,0)

! Initialize Fast Fourier Transform package
CALL SUFFT
CALL GSTATS(1802,1)



IF (LHOOK) CALL DR_HOOK('SETUP_TRANS',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------
9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

!endif INTERFACE

END SUBROUTINE SETUP_TRANS


