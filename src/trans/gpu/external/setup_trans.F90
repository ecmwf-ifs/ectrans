#define ALIGN(I, A) (((I)+(A)-1)/(A)*(A))
! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! (C) Copyright 2022- NVIDIA.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE SETUP_TRANS(KSMAX,KDGL,KDLON,KLOEN,LDSPLIT,PSTRET,&
 &                     KTMAX,KRESOL,PWEIGHT,LDGRIDONLY,LDUSERPNM,LDKEEPRPNM,LDUSEFLT,&
 &                     LDSPSETUPONLY,LDPNMONLY,LDUSEFFTW,&
 &                     LDLL,LDSHIFTLL,CDIO_LEGPOL,CDLEGPOLFNAME,KLEGPOLPTR,KLEGPOLPTR_LEN)

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

!     Explicit arguments : KLOEN,LDSPLIT are optional arguments
!     --------------------
!     KSMAX - spectral truncation required
!     KDGL  - number of Gaussian latitudes
!     KDLON - number of points on each Gaussian latitude [2*KDGL]
!     KLOEN(:) - number of points on each Gaussian latitude [2*KDGL]
!     LDSPLIT - true if split latitudes in grid-point space [false]
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
!     LDUSERPNM  - Use Belusov algorithm to compute legendre pol. (else new alg.)
!     LDKEEPRPNM - Keep Legendre Polynomials (only applicable when using
!                  FLT, otherwise always kept)
!     LDPNMONLY  - Compute the Legendre polynomials only, not the FFTs.
!     LDUSEFFTW    - Use FFTW for FFTs
!     LDLL                 - Setup second set of input/output latitudes
!                                 the number of input/output latitudes to transform is equal KDGL
!                                 or KDGL+2 in the case that includes poles + equator
!                                 the number of input/output longitudes to transform is 2*KDGL
!     LDSHIFTLL       - Shift output lon/lat data by 0.5*dx and 0.5*dy
!     CDIO_LEGPOL  - IO option on Legendre polinomials :  N.B. Only works for NPROC=1
!                    Options:
!                    'READF' -  read Leg.Pol. from file CDLEGPOLFNAME
!                    'WRITEF' - write Leg.Pol. to file CDLEGPOLFNAME
!                    'MEMBUF' - Leg. Pol provided in shared memory segment pointed to by KLEGPOLPTR of
!                               length KLEGPOLPTR_LEN
!     CDLEGPOLFNAME - file name for Leg.Pol. IO
!     KLEGPOLPTR    - pointer to Legendre polynomials memory segment
!     KLEGPOLPTR_LEN  - length of  Legendre polynomials memory segment

!     Method.
!     -------

!     Externals.  SET_RESOL   - set resolution
!     ----------  SETUP_DIMS  - setup distribution independent dimensions
!                 SUMP_TRANS_PRELEG - first part of setup of distr. environment
!                 SULEG - Compute Legandre polonomial and Gaussian
!                         Latitudes and Weights
!                 SUMP_TRANS - Second part of setup of distributed environment
!                 SUFFT - setup for FFT
!                 SHAREDMEM_CREATE - create memory buffer for Leg.pol.

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        Daan Degrauwe : Mar 2012 E'-zone dimensions
!        R. El Khatib 09-Aug-2012 %LAM in GEOM_TYPE
!        R. El Khatib 14-Jun-2013 PSTRET, LDPNMONLY, LENABLED
!        G. Mozdzynski : Oct 2014 Support f
!        N. Wedi       : Apr 2015 Support dual set of lat/lon
!        G. Mozdzynski : Jun 2015 Support alternative FFTs to FFTW
!        M.Hamrud/W.Deconinck : July 2015 IO options for Legenndre polynomials
!        R. El Khatib 07-Mar-2016 Better flexibility for Legendre polynomials computation in stretched mode
!     ------------------------------------------------------------------

USE PARKIND1,        ONLY: JPIM, JPRB, JPRD
USE PARKIND_ECTRANS, ONLY: JPRBT

!ifndef INTERFACE

USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_INT, C_ASSOCIATED, C_SIZE_T, C_SIZEOF
USE EC_ENV_MOD,                  ONLY: EC_GETENV
USE TPM_GEN,                     ONLY: NOUT, MSETUP0, NCUR_RESOL, NDEF_RESOL, &
     &                                 NMAX_RESOL, NPRINTLEV, LENABLED, NERR
USE TPM_DIM,                     ONLY: R, DIM_RESOL
USE TPM_DISTR,                   ONLY: D, DISTR_RESOL, NPROC, NPRTRV, MYPROC
USE TPM_GEOMETRY,                ONLY: G, GEOM_RESOL
USE TPM_FIELDS,                  ONLY: FIELDS_RESOL, F
USE TPM_FIELDS_GPU,              ONLY: FIELDS_GPU_RESOL, FG
USE TPM_FLT,                     ONLY: FLT_RESOL, S
USE TPM_CTL,                     ONLY: CTL_RESOL, C
USE SET_RESOL_MOD,               ONLY: SET_RESOL
USE SETUP_DIMS_MOD,              ONLY: SETUP_DIMS
USE SUMP_TRANS_MOD,              ONLY: SUMP_TRANS
USE SUMP_TRANS_PRELEG_MOD,       ONLY: SUMP_TRANS_PRELEG
USE SULEG_MOD,                   ONLY: SULEG
USE PRE_SULEG_MOD,               ONLY: PRE_SULEG
USE SUFFT_MOD,                   ONLY: SUFFT
USE ABORT_TRANS_MOD,             ONLY: ABORT_TRANS
USE SHAREDMEM_MOD,               ONLY: SHAREDMEM_CREATE
USE YOMHOOK,                     ONLY: LHOOK, DR_HOOK, JPHOOK
USE PREPSNM_MOD,                 ONLY: PREPSNM
#ifdef ACCGPU
USE OPENACC,                     ONLY: ACC_DEVICE_KIND, ACC_GET_DEVICE_TYPE, ACC_GET_NUM_DEVICES, &
  &                                    ACC_SET_DEVICE_NUM, ACC_GET_DEVICE_NUM
#endif
#ifdef OMPGPU
! TODO: add OMP equivalents to ACC library routines
!USE OMP_LIB
#endif

!endif INTERFACE

IMPLICIT NONE

! Dummy arguments

INTEGER(KIND=JPIM) ,INTENT(IN) :: KSMAX,KDGL
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN) :: KDLON
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN) :: KLOEN(:)
LOGICAL            ,OPTIONAL,INTENT(IN) :: LDSPLIT
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN) :: KTMAX
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(OUT):: KRESOL
REAL(KIND=JPRD)    ,OPTIONAL,INTENT(IN) :: PWEIGHT(:)
REAL(KIND=JPRD)    ,OPTIONAL,INTENT(IN) :: PSTRET
LOGICAL   ,OPTIONAL,INTENT(IN):: LDGRIDONLY
LOGICAL   ,OPTIONAL,INTENT(IN):: LDUSEFLT
LOGICAL   ,OPTIONAL,INTENT(IN):: LDUSERPNM
LOGICAL   ,OPTIONAL,INTENT(IN):: LDKEEPRPNM
LOGICAL   ,OPTIONAL,INTENT(IN):: LDSPSETUPONLY
LOGICAL   ,OPTIONAL,INTENT(IN):: LDPNMONLY
LOGICAL   ,OPTIONAL,INTENT(IN):: LDUSEFFTW
LOGICAL   ,OPTIONAL,INTENT(IN):: LDLL
LOGICAL   ,OPTIONAL,INTENT(IN):: LDSHIFTLL
CHARACTER(LEN=*),OPTIONAL,INTENT(IN):: CDIO_LEGPOL
CHARACTER(LEN=*),OPTIONAL,INTENT(IN):: CDLEGPOLFNAME
TYPE(C_PTR) ,OPTIONAL,INTENT(IN) :: KLEGPOLPTR
INTEGER(C_SIZE_T) ,OPTIONAL,INTENT(IN) :: KLEGPOLPTR_LEN

!ifndef INTERFACE

! Local variables
INTEGER(KIND=JPIM) :: JGL,JRES,IDEF_RESOL
INTEGER(KIND=JPIM) :: JMLOC, KM, ILA, ILS, KMLOC, KDGLU, JK, I, J

INTEGER(KIND=JPIM) :: IPROC, IPROCS, ISTAN, ISTAS, ISL, IGLS, JFLD, IMLOC0(1)

LOGICAL :: LLP1,LLP2, LLSPSETUPONLY
REAL(KIND=JPRD)    :: ZTIME0,ZTIME1,ZTIME2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

CHARACTER(LEN=8)  :: CENV

#ifdef ACCGPU
INTEGER(ACC_DEVICE_KIND) :: IDEVTYPE
#endif
INTEGER :: INUMDEVS, IUNIT, ISTAT, IDEV, MYGPU

#include "user_clock.intfb.h"
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SETUP_TRANS',0,ZHOOK_HANDLE)

IF(MSETUP0 == 0) THEN
  CALL ABORT_TRANS('SETUP_TRANS: SETUP_TRANS0 HAS TO BE CALLED BEFORE SETUP_TRANS')
ENDIF
LLP1 = NPRINTLEV>0
LLP2 = NPRINTLEV>1
IF(LLP1) WRITE(NOUT,*) '=== ENTER ROUTINE SETUP_TRANS ==='
IF(LLP1) THEN
  IF (JPRBT == JPRD) THEN
    WRITE(NOUT,'(A)') "GPU double precision version, with following compile-time options : "
  ELSE
    WRITE(NOUT,'(A)') "GPU single precision version, with following compile-time options : "
  ENDIF
#ifdef ACCGPU
  WRITE(NOUT,'(A)') " - OpenACC-based offload"
#else
  WRITE(NOUT,'(A)') " - OpenMP-based offload"
#endif
#ifdef USE_GPU_AWARE_MPI
  WRITE(NOUT,'(A)') " - GPU-aware MPI"
#endif
#ifdef USE_GRAPHS_FFT
  WRITE(NOUT,'(A)') " - graph-based FFT scheduling"
#endif
#ifdef USE_GRAPHS_GEMM
  WRITE(NOUT,'(A)') " - graph-based GEMM scheduling"
#endif
#ifdef USE_CUTLASS
  WRITE(NOUT,'(A)') " - Cutlass-based GEMM operations"
#endif
#ifdef USE_3XTF32
  WRITE(NOUT,'(A)') " - tensor-core usage for 32b Cutlass operations"
#endif
  WRITE(NOUT,'(A)')
ENDIF


! Allocate resolution dependent structures
IF(.NOT. ALLOCATED(DIM_RESOL)) THEN
  IDEF_RESOL = 1
  ALLOCATE(DIM_RESOL(NMAX_RESOL))
  ALLOCATE(FIELDS_RESOL(NMAX_RESOL))
  ALLOCATE(FIELDS_GPU_RESOL(NMAX_RESOL))
  ALLOCATE(GEOM_RESOL(NMAX_RESOL))
  ALLOCATE(DISTR_RESOL(NMAX_RESOL))
  ALLOCATE(FLT_RESOL(NMAX_RESOL))
  ALLOCATE(CTL_RESOL(NMAX_RESOL))
  GEOM_RESOL(:)%LAM=.FALSE.
  ALLOCATE(LENABLED(NMAX_RESOL))
  LENABLED(:)=.FALSE.
ELSE
  IDEF_RESOL = NMAX_RESOL+1
  DO JRES=1,NMAX_RESOL
    IF(.NOT.LENABLED(JRES)) THEN
      IDEF_RESOL = JRES
      EXIT
    ENDIF
  ENDDO
  IF(IDEF_RESOL > NMAX_RESOL) THEN
    CALL ABORT_TRANS('SETUP_TRANS:IDEF_RESOL > NMAX_RESOL')
  ENDIF
ENDIF

IF (PRESENT(KRESOL)) THEN
  KRESOL=IDEF_RESOL
ENDIF

! Point at structures due to be initialized
CALL SET_RESOL(IDEF_RESOL,LDSETUP=.TRUE.)

IF(LLP1) WRITE(NOUT,*) '=== DEFINING RESOLUTION ',NCUR_RESOL



! Defaults for optional arguments


G%LREDUCED_GRID = .FALSE.
G%RSTRET=1.0_JPRBT
D%LGRIDONLY = .FALSE.
D%LSPLIT = .FALSE.
D%LCPNMONLY=.FALSE.
S%LUSE_BELUSOV=.TRUE. ! use Belusov algorithm to compute RPNM array instead of per m
S%LKEEPRPNM=.FALSE. ! Keep Legendre polonomials (RPNM)
LLSPSETUPONLY = .FALSE. ! Only create distributed spectral setup
S%LDLL = .FALSE. ! use mapping to/from second set of latitudes
S%LSHIFTLL = .FALSE. ! shift output lat-lon by 0.5dx, 0.5dy
C%LREAD_LEGPOL = .FALSE.
C%LWRITE_LEGPOL = .FALSE.


! NON-OPTIONAL ARGUMENTS
R%NSMAX = KSMAX
R%NDGL  = KDGL
! E'-defaults
R%NNOEXTZL=0
R%NNOEXTZG=0


IF(PRESENT(LDSPSETUPONLY)) THEN
  LLSPSETUPONLY=LDSPSETUPONLY
! <<<<<<<<<<< EXTRA TO WORKAROUND NOT YET IMPLEMENTED FEATURE
  IF (LLSPSETUPONLY) THEN
    WRITE(NOUT,'(A)') "DEVELOPER WARNING: LDSPSETUPONLY IS NOT YET IMPLEMENTED CORRECTLY WITH GPU BACKEND. IGNORING IT FOR NOW"
    LLSPSETUPONLY = .FALSE.
    R%NDGL = NPROC
    ! Make even and positive
    IF (MOD(R%NDGL,2) /= 0) THEN
      R%NDGL = NPROC+1
    ENDIF
    R%NDGL = MAX(2,R%NDGL)
  ENDIF
! >>>>>>>>>>>>>
ENDIF


! IMPLICIT argument :
G%LAM = .FALSE.

IF(PRESENT(KDLON)) THEN
  R%NDLON = KDLON
ELSE
  R%NDLON = 2*R%NDGL
ENDIF

IF(PRESENT(LDLL)) THEN
  S%LDLL=LDLL
  IF( LDLL ) THEN
    CALL ABORT_TRANS ('SETUP_TRANS: LDLL=.TRUE. is not yet supported with GPU backend')

    S%NDLON=R%NDLON
    ! account for pole + equator
    R%NDGL=R%NDGL+2
    IF(PRESENT(LDSHIFTLL)) THEN
      S%LSHIFTLL = LDSHIFTLL
      ! geophysical (shifted) lat-lon without pole and equator
      IF(S%LSHIFTLL) R%NDGL=R%NDGL-2
    ENDIF
    S%NDGL=R%NDGL
  ENDIF
ENDIF

IF (R%NDGL <= 0 .OR. MOD(R%NDGL,2) /= 0) THEN
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

IF(PRESENT(PWEIGHT)) THEN
  D%LWEIGHTED_DISTR = .TRUE.
  IF( D%LWEIGHTED_DISTR .AND. .NOT.D%LSPLIT )THEN
    CALL ABORT_TRANS('SETUP_TRANS: LWEIGHTED_DISTR=T AND LSPLIT=F NOT SUPPORTED')
  ENDIF
  IF(SIZE(PWEIGHT) /= SUM(G%NLOEN(:)) )THEN
    CALL ABORT_TRANS('SETUP_TRANS:SIZE(PWEIGHT) /= SUM(G%NLOEN(:))')
  ENDIF
  IF( MINVAL(PWEIGHT(:)) < 0.0_JPRBT )THEN
    CALL ABORT_TRANS('SETUP_TRANS: INVALID WEIGHTS')
  ENDIF
  ALLOCATE(D%RWEIGHT(SIZE(PWEIGHT)))
  D%RWEIGHT(:)=PWEIGHT(:)
ELSE
  D%LWEIGHTED_DISTR = .FALSE.
ENDIF

IF(PRESENT(LDGRIDONLY)) THEN
  D%LGRIDONLY=LDGRIDONLY
! <<<<<<<<<<< EXTRA TO WORKAROUND NOT YET IMPLEMENTED FEATURE
  IF (D%LGRIDONLY) THEN
      R%NSMAX=1
      R%NTMAX = R%NSMAX
      WRITE(NOUT,'(A,I0)') "DEVELOPER WARNING: LDGRIDONLY IS NOT YET IMPLEMENTED CORRECTLY WITH GPU BACKEND. IGNORE AND USE TRUNCATION: ", R%NSMAX
      D%LGRIDONLY = .FALSE.
  ENDIF
! >>>>>>>>>>>>>
ENDIF

IF(PRESENT(LDPNMONLY)) THEN
  D%LCPNMONLY=LDPNMONLY
ENDIF

! Setup distribution independent dimensions
CALL SETUP_DIMS

S%LSOUTHPNM=.FALSE.
IF(PRESENT(PSTRET)) THEN
  IF (ABS(PSTRET-1.0_JPRBT)>100._JPRBT*EPSILON(1._JPRBT)) THEN
    G%RSTRET=PSTRET
    S%LSOUTHPNM=.TRUE.
    R%NLEI3=2*R%NLEI3 ! double
  ENDIF
ENDIF

IF(PRESENT(CDIO_LEGPOL)) THEN
  IF(NPROC > 1) CALL  ABORT_TRANS('SETUP_TRANS:CDIO_LEGPOL OPTIONS ONLY FOR NPROC=1 ')
  IF(TRIM(CDIO_LEGPOL) == 'readf' .OR. TRIM(CDIO_LEGPOL) == 'READF' ) THEN
    IF(.NOT.PRESENT(CDLEGPOLFNAME)) CALL  ABORT_TRANS('SETUP_TRANS: CDLEGPOLFNAME ARGUMENT MISSING')
    C%LREAD_LEGPOL = .TRUE.
    C%CLEGPOLFNAME = TRIM(CDLEGPOLFNAME)
    C%CIO_TYPE='file'
  ELSEIF(TRIM(CDIO_LEGPOL) == 'writef' .OR. TRIM(CDIO_LEGPOL) == 'WRITEF') THEN
    IF(.NOT.PRESENT(CDLEGPOLFNAME)) CALL  ABORT_TRANS('SETUP_TRANS: CDLEGPOLFNAME ARGUMENT MISSING')
    C%LWRITE_LEGPOL = .TRUE.
    C%CLEGPOLFNAME = TRIM(CDLEGPOLFNAME)
    C%CIO_TYPE='file'
  ELSEIF(TRIM(CDIO_LEGPOL) == 'membuf' .OR. TRIM(CDIO_LEGPOL) == 'MEMBUF') THEN
    IF(.NOT.PRESENT(KLEGPOLPTR)) CALL  ABORT_TRANS('SETUP_TRANS: KLEGPOLPTR  ARGUMENT MISSING')
    IF(.NOT.C_ASSOCIATED(KLEGPOLPTR))  CALL  ABORT_TRANS('SETUP_TRANS: KLEGPOLPTR NULL POINTER')
    IF(.NOT.PRESENT(KLEGPOLPTR_LEN)) CALL  ABORT_TRANS('SETUP_TRANS: KLEGPOLPTR_LEN ARGUMENT MISSING')
    C%LREAD_LEGPOL = .TRUE.
    C%CIO_TYPE='mbuf'
    CALL SHAREDMEM_CREATE( C%STORAGE,KLEGPOLPTR,KLEGPOLPTR_LEN)
  ELSE
    WRITE(NERR,*) 'CDIO_LEGPOL ', TRIM(CDIO_LEGPOL)
    CALL  ABORT_TRANS('SETUP_TRANS:CDIO_LEGPOL UNKNOWN METHOD ')
  ENDIF
ENDIF

IF(PRESENT(LDUSEFLT)) THEN
  IF (LDUSEFLT) THEN
    CALL ABORT_TRANS('SETUP_TRANS: LDUSEFLT option is not supported for GPU')
  ENDIF
ENDIF
IF(PRESENT(LDUSERPNM)) THEN
  S%LUSE_BELUSOV=LDUSERPNM
ENDIF
IF(PRESENT(LDKEEPRPNM)) THEN
  S%LKEEPRPNM=LDKEEPRPNM
ENDIF
!     Setup resolution dependent structures
!     -------------------------------------

! First part of setup of distributed environment
CALL SUMP_TRANS_PRELEG

IF( .NOT.LLSPSETUPONLY ) THEN

! Compute Legendre polonomial and Gaussian Latitudes and Weights
  CALL SULEG

! Second part of setup of distributed environment
  CALL SUMP_TRANS
  CALL GSTATS(1802,0)

! Initialize Fast Fourier Transform package
  IF (.NOT.D%LCPNMONLY) CALL SUFFT
  CALL GSTATS(1802,1)
ELSE
  CALL PRE_SULEG
ENDIF

! Signal the current resolution is active
LENABLED(IDEF_RESOL)=.TRUE.
NDEF_RESOL = COUNT(LENABLED)

IF (LHOOK) CALL DR_HOOK('SETUP_TRANS',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------
9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

IF( .NOT.D%LGRIDONLY ) THEN

IUNIT=300+MYPROC

#ifdef ACCGPU
!!IDEVTYPE=ACC_DEVICE_NVIDIA
IDEVTYPE=ACC_GET_DEVICE_TYPE()
INUMDEVS = ACC_GET_NUM_DEVICES(IDEVTYPE)
MYGPU = MOD(MYPROC-1,INUMDEVS)
CALL ACC_SET_DEVICE_NUM(MYGPU, IDEVTYPE)
MYGPU = ACC_GET_DEVICE_NUM(IDEVTYPE)
!ISTAT  = CUDA_GETDEVICE(IDEV)
#endif

WRITE(NOUT,*) 'R%NTMAX=',R%NTMAX
WRITE(NOUT,*) 'R%NSMAX=',R%NSMAX

! Initialize A arrays

ALLOCATE(FG%ZAA(ALIGN(R%NDGNH,8),ALIGN((R%NTMAX+2)/2,8),D%NUMP))
ALLOCATE(FG%ZAS(ALIGN(R%NDGNH,8),ALIGN((R%NTMAX+3)/2,8),D%NUMP))

FG%ZAA(:,:,:) = 0._JPRBT
FG%ZAS(:,:,:) = 0._JPRBT

DO JMLOC=1,D%NUMP
  KM = D%MYMS(JMLOC)
  KDGLU = G%NDGLU(KM)
  ILA = (R%NSMAX-KM+2)/2
  ILS = (R%NSMAX-KM+3)/2

  FG%ZAA(1:KDGLU,1:ILA,JMLOC)=S%FA(JMLOC)%RPNMA(1:KDGLU,1:ILA)
  FG%ZAS(1:KDGLU,1:ILS,JMLOC)=S%FA(JMLOC)%RPNMS(1:KDGLU,1:ILS)
ENDDO

! arrays for m=0 in ledir_mod:
IMLOC0 = FINDLOC(D%MYMS,0)
IF(IMLOC0(1) > 0) THEN
  ALLOCATE(FG%ZAA0(SIZE(FG%ZAA,1),SIZE(FG%ZAA,2)))
  ALLOCATE(FG%ZAS0(SIZE(FG%ZAS,1),SIZE(FG%ZAS,2)))
  FG%ZAA0 = FG%ZAA(:,:,IMLOC0(1))
  FG%ZAS0 = FG%ZAS(:,:,IMLOC0(1))
ENDIF

ALLOCATE(FG%ZEPSNM(D%NUMP,0:R%NTMAX+2))
FG%ZEPSNM = 0._JPRBT
CALL PREPSNM !Initialize on the host

WRITE(NOUT,*)'setup_trans: sizes1 NUMP=',D%NUMP
#ifdef ACCGPU
  WRITE(NOUT,*) 'Using OpenACC'
#endif
#ifdef OMPGPU
  WRITE(NOUT,*) 'Using OpenMP offloading'
#endif
WRITE(NOUT,'(A10,":",I9,"B")') 'FG%ZAS', C_SIZEOF(FG%ZAS(1,1,1))*SIZE(FG%ZAS)
WRITE(NOUT,'(A10,":",I9,"B")') 'FG%ZAA', C_SIZEOF(FG%ZAA(1,1,1))*SIZE(FG%ZAA)
WRITE(NOUT,'(A10,":",I9,"B")') 'FG%ZAS0', C_SIZEOF(FG%ZAS0(1,1))*SIZE(FG%ZAS0)
WRITE(NOUT,'(A10,":",I9,"B")') 'FG%ZAA0', C_SIZEOF(FG%ZAA0(1,1))*SIZE(FG%ZAA0)
WRITE(NOUT,'(A10,":",I9,"B")') 'FG%ZEPSNM', C_SIZEOF(FG%ZEPSNM(1,1))*SIZE(FG%ZEPSNM)

IF (IMLOC0(1) > 0) THEN
#ifdef ACCGPU
  !$ACC ENTER DATA COPYIN(FG%ZAA0,FG%ZAS0) ASYNC(1)
#endif
#ifdef OMPGPU
  !$OMP TARGET ENTER DATA MAP(TO:FG%ZAA0,FG%ZAS0)
#endif
ENDIF
#ifdef ACCGPU
#ifdef _CRAYFTN
!$ACC ENTER DATA COPYIN(R,R%NSMAX,R%NTMAX,R%NDGL,R%NDGNH) ASYNC(1)
#else
!$ACC ENTER DATA COPYIN(R) ASYNC(1)
#endif
!$ACC ENTER DATA COPYIN(F,F%RLAPIN,F%RACTHE,F%RW) ASYNC(1)
!$ACC ENTER DATA COPYIN(FG,FG%ZAA,FG%ZAS,FG%ZEPSNM) ASYNC(1)
#ifdef _CRAYFTN
!$ACC ENTER DATA COPYIN(D,D%NUMP,D%MYMS,D%NPNTGTB0,D%NPNTGTB1,D%NSTAGT0B,D%NSTAGT1B,D%NSTAGTF,D%NPROCM,D%NPROCL)&
!$ACC&           COPYIN(D%NPTRLS,D%MSTABF,D%NASM0,D%OFFSETS_GEMM1,D%OFFSETS_GEMM2,D%NDGL_FS) ASYNC(1)
#else
!$ACC ENTER DATA COPYIN(D,D%MYMS,D%NPNTGTB0,D%NPNTGTB1,D%NSTAGT0B,D%NSTAGT1B,D%NSTAGTF,D%NPROCM,D%NPROCL)&
!$ACC&           COPYIN(D%NPTRLS,D%MSTABF,D%NASM0,D%OFFSETS_GEMM1,D%OFFSETS_GEMM2) ASYNC(1)
#endif
!$ACC ENTER DATA COPYIN(G,G%NDGLU,G%NMEN,G%NLOEN) ASYNC(1)
!$ACC WAIT(1)
#endif
#ifdef OMPGPU
!$OMP TARGET ENTER DATA MAP(ALLOC:FG%ZAA,FG%ZAS)
!$OMP TARGET ENTER DATA MAP(TO:FG,F,S,D,R,G)
!$OMP BARRIER
#endif

WRITE(NOUT,*) '===GPU arrays successfully allocated'

! TODO: This might be good idea - those polynomials are not needed
!DO JMLOC=1,D%NUMP
!  DEALLOCATE(S%FA(JMLOC)%RPNMA)
!  DEALLOCATE(S%FA(JMLOC)%RPNMS)
!ENDDO

!endif INTERFACE

ENDIF ! D%LGRIDONLY

END SUBROUTINE SETUP_TRANS
