! (C) Copyright 2000- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE SETUP_TRANS(KSMAX,KDGL,KDLON,KLOEN,LDSPLIT,PSTRET,&
&KFLEV,KTMAX,KRESOL,PWEIGHT,LDGRIDONLY,LDUSERPNM,LDKEEPRPNM,LDUSEFLT,&
&LDSPSETUPONLY,LDPNMONLY,LDUSEFFTW,&
&LDLL,LDSHIFTLL,CDIO_LEGPOL,CDLEGPOLFNAME,KLEGPOLPTR,KLEGPOLPTR_LEN)

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

USE PARKIND1        ,ONLY : JPIM     ,JPRB ,  JPRD
USE PARKIND_ECTRANS ,ONLY : JPRBT
USE, INTRINSIC :: ISO_C_BINDING, ONLY:  C_PTR, C_INT,C_ASSOCIATED,C_SIZE_T

!ifndef INTERFACE

USE TPM_GEN         ,ONLY : NOUT, MSETUP0, NCUR_RESOL, NDEF_RESOL, &
     &                      NMAX_RESOL, NPRINTLEV, LENABLED, NERR
USE TPM_DIM         ,ONLY : R, DIM_RESOL, R_NSMAX,R_NTMAX, R_NDGNH, R_NDGL
USE TPM_DISTR       ,ONLY : D, DISTR_RESOL,NPROC,nprtrv, D_NUMP,D_MYMS,D_NSTAGT0B,D_NSTAGT1B,D_NPROCL,D_NPNTGTB1, D_NASM0, &
& D_NSTAGTF,D_MSTABF,D_NPNTGTB0,D_NPROCM,D_NPTRLS,mysetv,mysetw, MYPROC
USE TPM_GEOMETRY    ,ONLY : G, GEOM_RESOL, G_NDGLU, G_NMEN, G_NMEN_MAX,G_NLOEN, G_NLOEN_MAX
USE TPM_FIELDS      ,ONLY : FIELDS_RESOL, F,F_RW, ZEPSNM, &
& ZAA,ZAS,TDZAA,TDZAS,&
& ZAA0,&
& ZAS0,KMLOC0
USE TPM_FFT         ,ONLY : T, FFT_RESOL
#ifdef WITH_FFTW
USE TPM_FFTW        ,ONLY : TW, FFTW_RESOL
#endif
USE TPM_FLT
USE TPM_TRANS       ,ONLY : ZAVE, ZMINGL, ZMAXGL, ZMINGPN, ZMAXGPN
USE TPM_CTL

USE SET_RESOL_MOD   ,ONLY : SET_RESOL
USE SETUP_DIMS_MOD  ,ONLY : SETUP_DIMS
USE SUMP_TRANS_MOD  ,ONLY : SUMP_TRANS
USE SUMP_TRANS_PRELEG_MOD ,ONLY : SUMP_TRANS_PRELEG
USE SULEG_MOD       ,ONLY : SULEG
USE PRE_SULEG_MOD   ,ONLY : PRE_SULEG
USE SUFFT_MOD       ,ONLY : SUFFT
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
USE SHAREDMEM_MOD    ,ONLY : SHAREDMEM_CREATE
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK,  JPHOOK
USE CUDA_DEVICE_MOD
USE PREPSNM_MOD     ,ONLY : PREPSNM
#ifdef _OPENACC
use openacc
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
REAL(KIND=JPRB)    ,OPTIONAL,INTENT(IN) :: PWEIGHT(:)
REAL(KIND=JPRB)    ,OPTIONAL,INTENT(IN) :: PSTRET
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN) :: KFLEV
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
INTEGER(KIND=JPIM),PARAMETER :: IMAXFLD=240
INTEGER(KIND=JPIM) :: JGL,JRES,IDEF_RESOL
INTEGER(KIND=JPIM) :: JMLOC, KM, ILA, ILS, KMLOC, KDGLU, JK, i, J

LOGICAL :: LLP1,LLP2, LLSPSETUPONLY
REAL(KIND=JPRD)    :: ZTIME0,ZTIME1,ZTIME2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

integer :: idevtype, inumdevs, mygpu, iunit, istat, idev

#include "user_clock.intfb.h"
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
  IDEF_RESOL = 1
  ALLOCATE(DIM_RESOL(NMAX_RESOL))
  ALLOCATE(FIELDS_RESOL(NMAX_RESOL))
  ALLOCATE(GEOM_RESOL(NMAX_RESOL))
  ALLOCATE(DISTR_RESOL(NMAX_RESOL))
  ALLOCATE(FFT_RESOL(NMAX_RESOL))
#ifdef WITH_FFTW
  ALLOCATE(FFTW_RESOL(NMAX_RESOL))
#endif
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
S%LUSEFLT=.FALSE. ! Use fast legendre transforms
#ifdef WITH_FFTW
TW%LFFTW=.FALSE. ! Use FFTW interface for FFTs
#endif
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
ENDIF

IF(PRESENT(LDSPSETUPONLY)) THEN
  LLSPSETUPONLY=LDSPSETUPONLY
ENDIF

IF(PRESENT(LDPNMONLY)) THEN
  D%LCPNMONLY=LDPNMONLY
ENDIF


#ifdef WITH_FFTW
IF(PRESENT(LDUSEFFTW)) THEN
  TW%LFFTW=LDUSEFFTW
ENDIF
IF( LLSPSETUPONLY .OR. D%LGRIDONLY ) THEN
  TW%LFFTW = .FALSE.
ENDIF
#endif

S%LSOUTHPNM=.FALSE.
IF(PRESENT(PSTRET)) THEN
  IF (ABS(PSTRET-1.0_JPRBT)>100._JPRBT*EPSILON(1._JPRBT)) THEN
    G%RSTRET=PSTRET
    S%LSOUTHPNM=.TRUE.
  ENDIF
ENDIF

IF(PRESENT(CDIO_LEGPOL)) THEN
  IF(NPROC > 1) CALL  ABORT_TRANS('SETUP_TRANS:CDIO_LEGPOL OPTIONS ONLY FOR NPROC=1 ')
  IF(R%NSMAX > 511 ) S%LUSEFLT = .TRUE. !To save IO and memory
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
  S%LUSEFLT=LDUSEFLT
ENDIF
IF(PRESENT(LDUSERPNM)) THEN
  S%LUSE_BELUSOV=LDUSERPNM
ENDIF
IF(PRESENT(LDKEEPRPNM)) THEN
  IF(S%LUSEFLT) THEN
    IF(LDKEEPRPNM.AND..NOT.LDUSERPNM) THEN
      CALL ABORT_TRANS('SETUP_TRANS: LDKEEPRPNM=true with LDUSERPNM=false')
    ENDIF
  ENDIF
  S%LKEEPRPNM=LDKEEPRPNM
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

iunit=300+myproc

#ifdef _OPENACC
!!idevtype=acc_device_nvidia
idevtype=acc_get_device_type()
inumdevs = acc_get_num_devices(idevtype)
mygpu = mod(MYPROC-1,inumdevs)
CALL acc_set_device_num(mygpu, idevtype)
mygpu = acc_get_device_num(idevtype)
istat  = cuda_GetDevice(idev)
WRITE(iunit,*) '===now going to allocate GPU arrays on processor: ', myproc, ' device = ', mygpu, ' ',idev, ' of ', inumdevs
#endif

!leading and trailing dimensions of A for symmetric and antisymmetric cases
! (same for ltinv and ltdir)
TDZAA=(R%NTMAX+2)/2
TDZAS=(R%NTMAX+3)/2
print*,'R%NTMAX=',R%NTMAX
print*,'R%NSMAX=',R%NSMAX

!$ACC ENTER DATA &
!$ACC& COPYIN(F,F%RN,F%RLAPIN,D,D%MYMS,R,G,G%NDGLU) &
!$ACC& COPYIN(D%NPNTGTB0,D%NPNTGTB1,D%NSTAGT0B,D%NSTAGT1B,D%NSTAGTF,G%NMEN,D%NPROCM,D%NPTRLS,G,D%MSTABF)

! Initialize A arrays

ALLOCATE(ZAA(R%NDGNH,TDZAA,D%NUMP))
ALLOCATE(ZAS(R%NDGNH,TDZAS,D%NUMP))

write(nout,*)'setup_trans: sizes1 NUMP=',D%NUMP
write(nout,*)'ZAS:',size(ZAS)
write(nout,*)'ZAA:',size(ZAA)

ZAA(:,:,:) = 0
DO JMLOC=1,D%NUMP
  KM = D%MYMS(JMLOC)   
  KDGLU = MIN(R%NDGNH,G%NDGLU(KM))
   
  ILA = (R%NSMAX-KM+2)/2
  DO JK=1,KDGLU
    DO J=1,ILA
      ZAA(JK,J,JMLOC)=S%FA(JMLOC)%RPNMA(JK,J)
    ENDDO
  ENDDO
ENDDO

ZAS(:,:,:) = 0
DO JMLOC=1,D%NUMP
  KM = D%MYMS(JMLOC)
  KDGLU = MIN(R%NDGNH,G%NDGLU(KM))

  ILS = (R%NSMAX-KM+3)/2
  DO JK=1,KDGLU
    DO J=1,ILS
      ZAS(JK,J,JMLOC)=S%FA(JMLOC)%RPNMS(JK,J)
    ENDDO
  ENDDO
ENDDO
!$ACC ENTER DATA COPYIN(ZAA,ZAS)

ALLOCATE(ZEPSNM(d%nump,0:R%NTMAX+2))
write(nout,*)'ZEPSNM :',size(ZEPSNM)
ZEPSNM = 0._JPRBT
CALL PREPSNM
!$ACC ENTER DATA COPYIN(ZEPSNM)

! TODO: I guess tose might be needed again
! add arrays for GPNORM1
!ALLOCATE(ZAVE(IF_FS,R%NDGL))
!ALLOCATE(ZMINGL(IF_FS,R%NDGL))
!ALLOCATE(ZMAXGL(IF_FS,R%NDGL))
!ALLOCATE(ZMINGPN(IF_FS))
!ALLOCATE(ZMAXGPN(IF_FS))

!ZAVE = 0._JPRBT
!ZMINGL = 0._JPRBT
!ZMAXGL = 0._JPRBT
!ZMINGPN = 0._JPRBT
!ZMAXGPN = 0._JPRBT
!!$ACC ENTER DATA COPYIN(ZAVE,ZMINGL,ZMAXGL,ZMINGPN,ZMAXGPN)

!set up flat copies of constant data
R_NSMAX=R%NSMAX
R_NTMAX=R%NTMAX
R_NDGNH=R%NDGNH
R_NDGL=R%NDGL


ALLOCATE(D_NSTAGT0B(SIZE(D%NSTAGT0B)))
ALLOCATE(D_NSTAGT1B(SIZE(D%NSTAGT1B)))
ALLOCATE(D_NPNTGTB0(0:SIZE(D%NPNTGTB0,1)-1,SIZE(D%NPNTGTB0,2)))
ALLOCATE(D_NPNTGTB1(SIZE(D%NPNTGTB1,1),SIZE(D%NPNTGTB1,2)))
ALLOCATE(D_MYMS(SIZE(D%MYMS)))
ALLOCATE(D_NPROCL(SIZE(D%NPROCL)))
ALLOCATE(D_NASM0(0:SIZE(D%NASM0)-1))
ALLOCATE(D_NSTAGTF(SIZE(D%NSTAGTF)))
ALLOCATE(D_MSTABF(SIZE(D%MSTABF)))
ALLOCATE(D_NPROCM(0:SIZE(D%NPROCM)-1))
ALLOCATE(D_NPTRLS(SIZE(D%NPTRLS)))

ALLOCATE(G_NDGLU(0:SIZE(G%NDGLU)-1))
ALLOCATE(G_NMEN(SIZE(G%NMEN)))
ALLOCATE(G_NLOEN(SIZE(G%NLOEN)))

ALLOCATE(F_RW(SIZE(F%RW)))


DO I=0,SIZE(G%NDGLU)-1
   G_NDGLU(I)=G%NDGLU(I)
end DO

G_NMEN_MAX=0
DO I=1,SIZE(G%NMEN)
   G_NMEN(I)=G%NMEN(I)
   if (G_NMEN(I) .gt. G_NMEN_MAX) G_NMEN_MAX=G_NMEN(I)
end DO

G_NLOEN_MAX=0
DO I=1,SIZE(G%NLOEN)
   G_NLOEN(I)=G%NLOEN(I)
   if (G_NLOEN(I) .gt. G_NLOEN_MAX) G_NLOEN_MAX=G_NLOEN(I)
end DO

DO I=1,SIZE(D%NSTAGT0B)
   D_NSTAGT0B(I)=D%NSTAGT0B(I)
END DO

DO I=1,SIZE(D%NSTAGT1B)
   D_NSTAGT1B(I)=D%NSTAGT1B(I)
END DO

DO I=1,SIZE(D%NPROCL)
   D_NPROCL(I)=D%NPROCL(I)
END DO

DO I=0,SIZE(D%NASM0)-1
   D_NASM0(I)=D%NASM0(I)
END DO

DO I=1,SIZE(D%NSTAGTF)
   D_NSTAGTF(I)=D%NSTAGTF(I)
END DO

DO I=1,SIZE(D%MSTABF)
   D_MSTABF(I)=D%MSTABF(I)
END DO

DO I=0,SIZE(D%NPROCM)-1
   D_NPROCM(I)=D%NPROCM(I)
END DO

DO I=1,SIZE(D%NPTRLS)
   D_NPTRLS(I)=D%NPTRLS(I)
END DO

DO I=1,SIZE(D%NPNTGTB0,2)
   DO J=0,SIZE(D%NPNTGTB0,1)-1
      D_NPNTGTB0(J,I)=D%NPNTGTB0(J,I)
   end DO
END DO

DO I=1,SIZE(D%NPNTGTB1,2)
   DO J=1,SIZE(D%NPNTGTB1,1)
      D_NPNTGTB1(J,I)=D%NPNTGTB1(J,I)
   end DO
END DO

D_NUMP=D%NUMP

KMLOC0 = -1
DO I=1,SIZE(D%MYMS)
   D_MYMS(I)=D%MYMS(I)
   IF(D_MYMS(I) == 0) KMLOC0 = I
end DO

! arrays for m=0 in ledir_mod:
IF(KMLOC0 >= 0) THEN
  ALLOCATE(ZAA0(R%NDGNH,TDZAA))
  ALLOCATE(ZAS0(R%NDGNH,TDZAS))
  ZAA0 = ZAA(:,:,KMLOC0)
  ZAS0 = ZAS(:,:,KMLOC0)
  !$ACC ENTER DATA COPYIN(ZAA0,ZAS0)
  WRITE(NOUT,*) 'GPU arrays for m=0 successfully allocated'
ENDIF

DO I=1,SIZE(F%RW)
   F_RW(I)=F%RW(I)
END DO

!$ACC ENTER DATA COPYIN(D_NSTAGT0B,D_NSTAGT1B,&
!$ACC&                  D_NPNTGTB1,D_NPROCL,D_MYMS,D_NASM0,D_NSTAGTF,D_MSTABF,&
!$ACC&                  D_NPNTGTB0,D_NPROCM,D_NPTRLS,G_NDGLU,G_NMEN,G_NLOEN,&
!$ACC&                  F_RW)

WRITE(NOUT,*) '===GPU arrays successfully allocated'
!$ACC wait

! free memory
!DO JMLOC=1,D%NUMP
!  DEALLOCATE(S%FA(JMLOC)%RPNMA)
!  DEALLOCATE(S%FA(JMLOC)%RPNMS)
!ENDDO

!endif INTERFACE

ENDIF

END SUBROUTINE SETUP_TRANS
