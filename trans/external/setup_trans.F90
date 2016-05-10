SUBROUTINE SETUP_TRANS(KSMAX,KDGL,KLOEN,LDLINEAR_GRID,LDSPLIT,PSTRET,&
&KTMAX,KRESOL,PWEIGHT,LDGRIDONLY,LDUSERPNM,LDKEEPRPNM,LDUSEFLT,&
&LDSPSETUPONLY,LDPNMONLY,LDUSEFFT992,LDUSEFFTW)

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

!     Explicit arguments : KLOEN,LDLINEAR_GRID,LDSPLIT are optional arguments
!     --------------------
!     KSMAX - spectral truncation required
!     KDGL  - number of Gaussian latitudes
!     KLOEN(:) - number of points on each Gaussian latitude [2*KDGL]
!     LDSPLIT - true if split latitudes in grid-point space [false]
!     LDLINEAR_GRID - true if linear grid
!     KTMAX - truncation order for tendencies?
!     KRESOL - the resolution identifier
!     PWEIGHT - the weight per grid-point (for a weighted distribution)
!     LDGRIDONLY - true if only grid space is required

!     KSMAX,KDGL,KTMAX and KLOEN are GLOBAL variables desribing the resolution
!     in spectral and grid-point space

!     LDSPLIT describe the distribution among processors of grid-point data and
!     has no relevance if you are using a single processor
 
!     PSTRET     - stretching factor - for the case the Legendre polynomials are
!                  computed on the stretched sphere - works with LSOUTHPNM
!     LDUSEFLT   - use Fast Legandre Transform (Butterfly algorithm)
!     LDUSERPNM  - Use Belusov to compute legendre pol. (else new alg.)
!     LDKEEPRPNM - Keep Legendre Polynomials (only applicable when using
!                  FLT, otherwise always kept)
!     LDPNMONLY  - Compute the Legendre polynomials only, not the FFTs.
!     LDUSEFFT992 - Use FF992 for FFTs (default)
!     LDUSEFFTW   - Use FFTW for FFTs
 
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
!        Daan Degrauwe : Mar 2012 E'-zone dimensions
!        R. El Khatib 09-Aug-2012 %LAM in GEOM_TYPE
!        R. El Khatib 14-Jun-2013 PSTRET, LDPNMONLY, LENABLED
!        G. Mozdzynski : Oct 2014 Support for FFTW
!        R. El Khatib 07-Mar-2016 Better flexibility for Legendre polynomials computation in stretched mode
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

!ifndef INTERFACE

USE TPM_GEN         ,ONLY : NOUT, MSETUP0, NCUR_RESOL, NDEF_RESOL, &
     &                      NMAX_RESOL, NPRINTLEV, LENABLED
USE TPM_DIM         ,ONLY : R, DIM_RESOL
USE TPM_DISTR       ,ONLY : D, DISTR_RESOL
USE TPM_GEOMETRY    ,ONLY : G, GEOM_RESOL
USE TPM_FIELDS      ,ONLY : FIELDS_RESOL
USE TPM_FFT         ,ONLY : T, FFT_RESOL
#ifdef WITH_FFTW
USE TPM_FFTW        ,ONLY : TW, FFTW_RESOL
#endif
USE TPM_FLT

USE SET_RESOL_MOD   ,ONLY : SET_RESOL
USE SETUP_DIMS_MOD  ,ONLY : SETUP_DIMS
USE SUMP_TRANS_MOD  ,ONLY : SUMP_TRANS
USE SUMP_TRANS_PRELEG_MOD ,ONLY : SUMP_TRANS_PRELEG
USE SULEG_MOD       ,ONLY : SULEG
USE SETUP_GEOM_MOD  ,ONLY : SETUP_GEOM
USE SUFFT_MOD       ,ONLY : SUFFT
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

!endif INTERFACE

IMPLICIT NONE

! Dummy arguments

INTEGER(KIND=JPIM) ,INTENT(IN) :: KSMAX,KDGL
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN) :: KLOEN(:)
LOGICAL   ,OPTIONAL,INTENT(IN) :: LDLINEAR_GRID
LOGICAL   ,OPTIONAL,INTENT(IN) :: LDSPLIT
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN) :: KTMAX
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(OUT):: KRESOL
REAL(KIND=JPRB)    ,OPTIONAL,INTENT(IN) :: PWEIGHT(:)
REAL(KIND=JPRB)    ,OPTIONAL,INTENT(IN) :: PSTRET
LOGICAL   ,OPTIONAL,INTENT(IN):: LDGRIDONLY
LOGICAL   ,OPTIONAL,INTENT(IN):: LDUSEFLT
LOGICAL   ,OPTIONAL,INTENT(IN):: LDUSERPNM
LOGICAL   ,OPTIONAL,INTENT(IN):: LDKEEPRPNM
LOGICAL   ,OPTIONAL,INTENT(IN):: LDSPSETUPONLY
LOGICAL   ,OPTIONAL,INTENT(IN):: LDPNMONLY
LOGICAL   ,OPTIONAL,INTENT(IN):: LDUSEFFT992
LOGICAL   ,OPTIONAL,INTENT(IN):: LDUSEFFTW

!ifndef INTERFACE

! Local variables
INTEGER(KIND=JPIM) :: JGL

LOGICAL :: LLP1,LLP2, LLSPSETUPONLY
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SETUP_TRANS',0,ZHOOK_HANDLE)

IF(MSETUP0 == 0) THEN
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
#ifdef WITH_FFTW
  ALLOCATE(FFTW_RESOL(NMAX_RESOL))
#endif
  ALLOCATE(FLT_RESOL(NMAX_RESOL))
  GEOM_RESOL(:)%LAM=.FALSE.
  ALLOCATE(LENABLED(NMAX_RESOL))
  LENABLED(:)=.FALSE.
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
G%RSTRET=1.0_JPRB
D%LGRIDONLY = .FALSE.
D%LSPLIT = .FALSE.
D%LCPNMONLY=.FALSE.
S%LUSERPNM=.TRUE. ! use RPNM array instead of per m
S%LKEEPRPNM=.TRUE. ! Keep Legendre polonomials (RPNM)
S%LUSEFLT=.FALSE. ! Use fast legendre transforms
T%LFFT992=.TRUE. ! Use FFT992 interface for FFTs
#ifdef WITH_FFTW
TW%LFFTW=.FALSE. ! Use FFTW interface for FFTs
#endif
LLSPSETUPONLY = .FALSE. ! Only create distributed spectral setup


! NON-OPTIONAL ARGUMENTS
R%NSMAX = KSMAX
R%NDGL  = KDGL
R%NDLON = 2*KDGL
! E'-defaults
R%NNOEXTZL=0
R%NNOEXTZG=0

! IMPLICIT argument :
G%LAM = .FALSE.

IF (KDGL <= 0 .OR. MOD(KDGL,2) /= 0) THEN
  CALL ABORT_TRANS ('SETUP_TRANS: KDGL IS NOT A POSITIVE, EVEN NUMBER')
ENDIF

! Optional arguments

ALLOCATE(G%NLOEN(R%NDGL))
IF(LLP2)WRITE(NOUT,9) 'NLOEN   ',SIZE(G%NLOEN   ),SHAPE(G%NLOEN   )
IF(PRESENT(KLOEN)) THEN
  IF( MINVAL(KLOEN(:)) <= 0 )THEN
     CALL ABORT_TRANS ('SETUP_TRANS: KLOEN INVALID (ONE or MORE POINTS <= 0)')
  ENDIF
  R%NDLON=MAXVAL(KLOEN(:))
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

IF(PRESENT(KTMAX)) THEN
  R%NTMAX = KTMAX
ELSE
  R%NTMAX = R%NSMAX
ENDIF
!Temporary?
IF(PRESENT(LDLINEAR_GRID)) THEN
  G%LINEAR_GRID = LDLINEAR_GRID
ELSEIF(R%NSMAX > (R%NDLON+3)/3) THEN
  G%LINEAR_GRID = .TRUE.
ENDIF

IF(PRESENT(PWEIGHT)) THEN
  D%LWEIGHTED_DISTR = .TRUE.
  IF( D%LWEIGHTED_DISTR .AND. .NOT.D%LSPLIT )THEN
    CALL ABORT_TRANS('SETUP_TRANS: LWEIGHTED_DISTR=T AND LSPLIT=F NOT SUPPORTED')
  ENDIF
  IF(SIZE(PWEIGHT) /= SUM(G%NLOEN(:)) )THEN
    CALL ABORT_TRANS('SETUP_TRANS:SIZE(PWEIGHT) /= SUM(G%NLOEN(:))')
  ENDIF
  IF( MINVAL(PWEIGHT(:)) < 0.0_JPRB )THEN
    CALL ABORT_TRANS('SETUP_TRANS: INVALID WEIGHTS')
  ENDIF
  ALLOCATE(D%RWEIGHT(SIZE(PWEIGHT)))
  D%RWEIGHT(:)=PWEIGHT(:)
ELSE
  D%LWEIGHTED_DISTR = .FALSE.
ENDIF

IF(PRESENT(LDGRIDONLY)) THEN
  D%LGRIDONLY=LDGRIDONLY
ENDIF

IF(PRESENT(LDSPSETUPONLY)) THEN
  LLSPSETUPONLY=LDSPSETUPONLY
ENDIF

IF(PRESENT(LDPNMONLY)) THEN
  D%LCPNMONLY=LDPNMONLY
ENDIF

IF(PRESENT(LDUSEFLT)) THEN
  S%LUSEFLT=LDUSEFLT
ENDIF
IF(PRESENT(LDUSERPNM)) THEN
  S%LUSERPNM=LDUSERPNM
ENDIF
IF(PRESENT(LDKEEPRPNM)) THEN
  IF(S%LUSEFLT) THEN
    IF(LDKEEPRPNM.AND..NOT.LDUSERPNM) THEN
      CALL ABORT_TRANS('SETUP_TRANS: LDKEEPRPNM=true with LDUSERPNM=false')
    ENDIF
  ENDIF
  S%LKEEPRPNM=LDKEEPRPNM
ENDIF

! Allow both FFT992 and FFTW initialisation for now,
! could be useful for debugging purposes
IF(PRESENT(LDUSEFFT992)) THEN
  T%LFFT992=LDUSEFFT992
ENDIF
#ifdef WITH_FFTW
IF(PRESENT(LDUSEFFTW)) THEN
  TW%LFFTW=LDUSEFFTW
ENDIF
#endif

S%LSOUTHPNM=.FALSE.
IF(PRESENT(PSTRET)) THEN
  IF (ABS(PSTRET-1.0_JPRB)>100._JPRB*EPSILON(1._JPRB)) THEN
    G%RSTRET=PSTRET
    S%LSOUTHPNM=.TRUE.
  ENDIF
ENDIF

!     Setup resolution dependent structures
!     -------------------------------------

! Setup distribution independent dimensions
CALL SETUP_DIMS

! First part of setup of distributed environment
CALL SUMP_TRANS_PRELEG

IF( .NOT.LLSPSETUPONLY ) THEN

! Compute Legendre polonomial and Gaussian Latitudes and Weights
CALL SULEG

CALL GSTATS(1802,0)

! Compute arrays related to grid-point geometry
!CALL SETUP_GEOM mover to SULEG
CALL GSTATS(1802,1)

! Second part of setup of distributed environment
CALL SUMP_TRANS
CALL GSTATS(1802,0)

! Initialize Fast Fourier Transform package
IF (.NOT.D%LCPNMONLY) CALL SUFFT
CALL GSTATS(1802,1)

ENDIF

! Signal the current resolution is active
LENABLED(NDEF_RESOL)=.TRUE.

IF (LHOOK) CALL DR_HOOK('SETUP_TRANS',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------
9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

!endif INTERFACE

END SUBROUTINE SETUP_TRANS
