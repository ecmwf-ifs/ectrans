! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

PROGRAM TRANSFORM_TEST

!
! Spectral transform test
!
! This test performs spectral to real and real to spectral transforms repeated in
! timed loop.
!
!
! Author : George Mozdzynski
!

USE PARKIND1  ,ONLY  : JPIM     ,JPRB, JPRD
USE OML_MOD ,ONLY : OML_MAX_THREADS
USE MPL_MPIF
USE MPL_MODULE
USE YOMGSTATS, ONLY: JPMAXSTAT

IMPLICIT NONE

! Number of points in top/bottom latitudes
INTEGER(KIND=JPIM), PARAMETER :: MIN_OCTA_POINTS = 20

INTEGER(KIND=JPIM) :: ISTACK, GETSTACKUSAGE
REAL(KIND=JPRB), DIMENSION(1)  :: ZMAXERR(5), ZERR(5)
REAL(KIND=JPRB) :: ZMAXERRG

INTEGER(KIND=JPIM) :: NERR,NLIN,INSF,NSMAX,NDGL,NQ
INTEGER(KIND=JPIM) :: NOUT,NOUTDUMP,NSPEC2,NGPTOT,NGPTOTG,IFLD,IFLDS,ICODE,IOUTSF,JROC,JB
INTEGER(KIND=JPIM) :: IERR,NSPEC2G,IRET,NTYPE,I
INTEGER(KIND=JPIM) :: JA,IB,JPRTRV
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NLOEN(:),ITO(:),NPRCIDS(:)
INTEGER(KIND=JPIM) :: MYPROC,JJ
INTEGER   :: JSTEP
REAL(KIND=JPRD)    :: ZTINIT,ZTLOOP,TIMEF, ZTSTEPMAX, ZTSTEPMIN, ZTSTEPAVG, ZTSTEPMED
REAL(KIND=JPRD)    :: ZTSTEPMAX1, ZTSTEPMIN1, ZTSTEPAVG1, ZTSTEPMED1
REAL(KIND=JPRD)    :: ZTSTEPMAX2, ZTSTEPMIN2, ZTSTEPAVG2, ZTSTEPMED2
REAL(KIND=JPRD),ALLOCATABLE :: ZTSTEP(:), ZTSTEP1(:), ZTSTEP2(:)
REAL(KIND=JPRB),ALLOCATABLE :: ZNORMSP(:),ZNORMSP1(:),ZNORMDIV(:),ZNORMDIV1(:)
REAL(KIND=JPRB),ALLOCATABLE :: ZNORMVOR(:),ZNORMVOR1(:),ZNORMT(:),ZNORMT1(:)
REAL(KIND=JPRB),ALLOCATABLE :: ZNORM(:),ZNORM1(:)
REAL(KIND=JPRD) :: ZAVEAVE(0:JPMAXSTAT)

! GRID-POINT SPACE DATA STRUCTURES
REAL(KIND=JPRB), ALLOCATABLE :: ZWINDS (:,:,:,:) ! Multilevel fields at t and t-dt
REAL(KIND=JPRB), ALLOCATABLE, TARGET :: ZGMV   (:,:,:,:) ! Multilevel fields at t and t-dt
REAL(KIND=JPRB), ALLOCATABLE, TARGET :: ZGMVS  (:,:,:)   ! Single level fields at t and t-dt

! SPECTRAL SPACE DATA STRUCTURES
REAL(KIND=JPRB), ALLOCATABLE :: ZSPVORG(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSPDIVG(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSPSPG(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSPTG(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE, TARGET :: SP3D(:,:,:)
REAL(KIND=JPRB), POINTER :: ZVOR(:,:) => NULL()
REAL(KIND=JPRB), POINTER :: ZDIV(:,:) => NULL()
REAL(KIND=JPRB), POINTER :: ZT(:,:,:) => NULL()
REAL(KIND=JPRB), ALLOCATABLE :: ZSP(:,:)

LOGICAL :: LSTACK
LOGICAL :: LUSERPNM, LKEEPRPNM, LUSEFLT
LOGICAL :: LTRACE_STATS,LSTATS_OMP, LSTATS_COMMS, LSTATS_MPL
LOGICAL :: LSTATS,LBARRIER_STATS, LBARRIER_STATS2, LDETAILED_STATS
LOGICAL :: LSTATS_ALLOC, LSYNCSTATS, LSTATSCPU, LSTATS_MEM
LOGICAL :: LXML_STATS
LOGICAL :: LFFTW
INTEGER(KIND=JPIM) :: NSTATS_MEM, NTRACE_STATS, NPRNT_STATS
LOGICAL :: LMPOFF
INTEGER(KIND=JPIM) :: ITERS=100

! Whether to print verbose output or not
LOGICAL :: VERBOSE = .FALSE.

REAL(KIND=JPRB) :: ZRA=6371229._JPRB

INTEGER(KIND=JPIM) :: NMAX_RESOL
INTEGER(KIND=JPIM) :: NPROMATR
INTEGER(KIND=JPIM) :: NCOMBFLEN

INTEGER(KIND=JPIM) :: NPROC
INTEGER(KIND=JPIM) :: NTHREAD
INTEGER(KIND=JPIM) :: NPRGPNS
INTEGER(KIND=JPIM) :: NPRGPEW
INTEGER(KIND=JPIM) :: NPRTRV
INTEGER(KIND=JPIM) :: NPRTRW
INTEGER(KIND=JPIM) :: NSPECRESMIN
INTEGER(KIND=JPIM) :: MYSETV
INTEGER(KIND=JPIM) :: MYSETW
INTEGER(KIND=JPIM) :: MP_TYPE
INTEGER(KIND=JPIM) :: MBX_SIZE

INTEGER(KIND=JPIM), ALLOCATABLE :: NUMLL(:), IVSET(:)
INTEGER(KIND=JPIM) :: IVSETSC(1)

INTEGER(KIND=JPIM) :: NFLEVG, NFLEVL
! SUMPINI
INTEGER(KIND=JPIM) :: ISQR
LOGICAL :: LSYNC_TRANS
LOGICAL :: LEQ_REGIONS


INTEGER(KIND=JPIM) :: NPROMA
INTEGER(KIND=JPIM) :: NGPBLKS
! LOCALS
INTEGER(KIND=JPIM) :: IPRTRV
INTEGER(KIND=JPIM) :: IPRTRW
INTEGER(KIND=JPIM) :: IPRUSED, ILEVPP, IREST, ILEV, JLEV, ILASTLEV

LOGICAL :: LLINFO

INTEGER(KIND=JPIM) :: NDIMGMV ! Third dim. of GMV "(NPROMA,NFLEVG,NDIMGMV,NGPBLKS)"
INTEGER(KIND=JPIM) :: NDIMGMVS ! Second dim. GMVS "(NPROMA,NDIMGMVS,NGPBLKS)"

! For processing command line arguments
CHARACTER(LEN=32) :: ARG

!===================================================================================================

#include "setup_trans0.h"
#include "setup_trans.h"
#include "inv_trans.h"
#include "dir_trans.h"
#include "dist_spec.h"
#include "gath_grid.h"
#include "trans_inq.h"
#include "specnorm.h"
#include "abor1.intfb.h"
#include "gstats_setup.intfb.h"

!===================================================================================================
! Initialize parameters
!===================================================================================================

NERR = 0
NOUT = 6
! Unit number for file to dump 2D fields to
NOUTDUMP = 7
! Max number of resolutions
NMAX_RESOL=37

! NPROMA for trans lib
NPROMATR=0
! Size of comm buffer
NCOMBFLEN=1800000
! EQ REGIONS flag
LEQ_REGIONS=.TRUE.
! Message Passing switch
LMPOFF=.FALSE.
! Activate barrier sync
LSYNC_TRANS=.true.
! Number of procs
NPROC=0
! Grid-point decomp
NPRGPNS=0
NPRGPEW=0
! Spectral decomp
NPRTRW=0
NPRTRV=0
! Minimum spectral resolution
! Used for controlling NPRTRW
NSPECRESMIN=80
! Message passing type
MP_TYPE=2
! Mailbox size
MBX_SIZE=150000000
! GSTATS statistics
LSTATS=.TRUE.
LDETAILED_STATS=.FALSE.
LSTATS_OMP=.FALSE.
LSTATS_COMMS=.FALSE.
LSTATS_MPL=.FALSE.
LBARRIER_STATS=.FALSE.
LBARRIER_STATS2=.FALSE.
LSTATSCPU=.FALSE.
LSYNCSTATS=.FALSE.
LXML_STATS=.FALSE.
LTRACE_STATS=.FALSE.
NSTATS_MEM=0
LSTATS_MEM=.FALSE.
LSTATS_ALLOC=.FALSE.
NTRACE_STATS=0
NPRNT_STATS=1
LUSERPNM=.FALSE.
LKEEPRPNM=.FALSE.
! Use fast Legendre transforms
LUSEFLT=.FALSE.
! Output stack info
LSTACK=.FALSE.
! Use FFTW
LFFTW=.TRUE.

! Default number of vertical levels
NFLEVG=137
! Number of 3D grid-point fields in GMV
NDIMGMV=9
! Number of 2D grid-point fields in GMVS
! surface pressure, north south der, east-west der
NDIMGMVS=3
! Set defaults for options
NLIN    = 0
NDGL    = 0
NQ      = 2

! Number of iterations for transform test
ITERS=10

! Locals
ILASTLEV = 0

!===================================================================================================
! Read command-line arguments
!===================================================================================================

DO I = 1, COMMAND_ARGUMENT_COUNT()
  CALL GET_COMMAND_ARGUMENT(I, ARG)

  SELECT CASE(ARG)
    ! Verbose output
    CASE("-v", "--verbose")
      VERBOSE = .TRUE.
    CASE DEFAULT
      CALL ABOR1("Unrecognized command-line option: " // ARG)
  END SELECT
END DO

!===================================================================================================
! Message passing setup
! Participating processors limited by -P option
!===================================================================================================

CALL MPL_INIT()
!IF( LSTATS ) CALL GSTATS(0,0)
ZTINIT=TIMEF()

NPROC= MPL_NPROC()
MYPROC = MPL_MYRANK()
NTHREAD= OML_MAX_THREADS()

! ONLY OUTPUT TO STDOUT ON PE 1
IF( NPROC > 1 ) THEN
  IF( MYPROC /= 1 ) THEN
    OPEN(UNIT=NOUT, FILE='/dev/null')
  ENDIF
ENDIF

IF(LDETAILED_STATS)THEN
  LSTATS_OMP=.TRUE.
  LSTATS_COMMS=.TRUE.
  LSTATS_MPL=.TRUE.
  LSTATSCPU=.TRUE.
  NPRNT_STATS=NPROC
!  LSTATS_MEM=.TRUE.
!  LSTATS_ALLOC=.TRUE.
ENDIF

!===================================================================================================

ALLOCATE(NPRCIDS(NPROC))
DO JJ=1,NPROC
  NPRCIDS(JJ) = JJ
ENDDO

IF( NPROC <= 1 ) LMPOFF=.TRUE.
! COMPUTE NPRGPNS and NPRGPEW
! THIS VERSION SELECTS MOST SQUARE-LIKE DISTRIBUTION
! THESE WILL CHANGE IF LEQ_REGIONS=.TRUE.
IF( NPROC == 0 ) NPROC = 1
ISQR=INT(SQRT(REAL(NPROC,JPRB)))
DO JA=ISQR,NPROC
  IB=NPROC/JA
  IF( JA*IB == NPROC ) THEN
    NPRGPNS=MAX(JA,IB)
    NPRGPEW=MIN(JA,IB)
    EXIT
  ENDIF
ENDDO

! FROM SUMPINI, ALTHOUGH THIS
! SHOULD BE SPECIFIED IN NAMELIST
IF( NSPECRESMIN==0 ) NSPECRESMIN=NPROC

! COMPUTE NPRTRV AND NPRTRW
! IF NOT PROVIDED IN NAMELIST
IF( NPRTRV > 0 .OR. NPRTRW > 0 ) THEN
  IF( NPRTRV == 0 ) NPRTRV=NPROC/NPRTRW
  IF( NPRTRW == 0 ) NPRTRW=NPROC/NPRTRV
  IF( NPRTRW*NPRTRV /= NPROC ) CALL ABOR1('TRANSFORM_TEST:NPRTRW*NPRTRV /= NPROC')
  IF( NPRTRW > NSPECRESMIN ) CALL ABOR1('TRANSFORM_TEST:NPRTRW > NSPECRESMIN')
ELSE
  DO JPRTRV=4,NPROC
    NPRTRV=JPRTRV
    NPRTRW=NPROC/NPRTRV
    IF( NPRTRV*NPRTRW /= NPROC ) CYCLE
    IF( NPRTRV > NPRTRW ) EXIT
    IF( NPRTRW > NSPECRESMIN ) CYCLE
    IF( NPRTRW <= NSPECRESMIN/(2*OML_MAX_THREADS()) ) EXIT
  ENDDO
  ! GO FOR APPROX SQUARE PARTITION FOR BACKUP
  IF( NPRTRV*NPRTRW /= NPROC .OR. NPRTRW > NSPECRESMIN .OR. NPRTRV > NPRTRW ) THEN
    ISQR=INT(SQRT(REAL(NPROC,JPRB)))
    DO JA=ISQR,NPROC
      IB=NPROC/JA
      IF (JA*IB == NPROC) THEN
        NPRTRW=MAX(JA,IB)
        NPRTRV=MIN(JA,IB)
        IF (NPRTRW > NSPECRESMIN ) CALL ABOR1('TRANSFORM_TEST:NPRTRW &
                                           & (approx square value) > NSPECRESMIN')
        EXIT
      ENDIF
    ENDDO
  ENDIF
ENDIF

! Create communicators for MPI groups
IF (.NOT.LMPOFF) THEN
  CALL MPL_GROUPS_CREATE(NPRTRW,NPRTRV)
ENDIF

IF (LMPOFF) THEN
  MYSETW=(MYPROC-1)/NPRTRV+1
  MYSETV=MOD(MYPROC-1,NPRTRV)+1
ELSE
  CALL MPL_CART_COORDS(MYPROC,MYSETW,MYSETV)
  ! Just checking for now...
  IPRTRV=MOD(MYPROC-1,NPRTRV)+1
  IPRTRW=(MYPROC-1)/NPRTRV+1
  IF (IPRTRV/=MYSETV .OR. IPRTRW/=MYSETW) THEN
    CALL ABOR1('TRANSFORM_TEST:Inconsistency when computing MYSETW and MYSETV')
  ENDIF
ENDIF

IF (.NOT.LMPOFF) THEN
  LLINFO=.FALSE.
  IF (MYPROC == 1) LLINFO=.TRUE.
  CALL MPL_BUFFER_METHOD(KMP_TYPE=MP_TYPE,KMBX_SIZE=MBX_SIZE,KPROCIDS=NPRCIDS,LDINFO=LLINFO)
ENDIF

! Determine number of local levels for Fourier and Legendre calculations
! based on the values of NFLEVG and NPRTRV
ALLOCATE(NUMLL(NPRTRV+1))

! Calculate remainder
IPRUSED=MIN(NFLEVG+1,NPRTRV)
ILEVPP=NFLEVG/NPRTRV
IREST=NFLEVG-ILEVPP*NPRTRV
DO JROC=1,NPRTRV
  IF(JROC <= IREST) THEN
    NUMLL(JROC)=ILEVPP+1
  ELSE
    NUMLL(JROC)=ILEVPP
  ENDIF
ENDDO
NUMLL(IPRUSED+1:NPRTRV+1)=0

NFLEVL=NUMLL(MYSETV)

IVSETSC(1)=IPRUSED

IFLD=0
IFLDS=0

ICODE = 0

!===================================================================================================
! Set resolution parameters
!===================================================================================================

! Spectral truncation
NSMAX = 79
NDGL = 2 * (NSMAX + 1)

! Calculate number of points at each latitude for octahedral grid
ALLOCATE(NLOEN(NDGL))

DO I = 1, NDGL / 2
  NLOEN(I) = MIN_OCTA_POINTS + 4 * (I - 1)
  NLOEN(NDGL - I + 1) = NLOEN(I)
END DO

!===================================================================================================
! Call ecTrans setup routines
!===================================================================================================

CALL SETUP_TRANS0(KOUT=NOUT,KERR=NERR,KPRINTLEV=MERGE(2, 0, VERBOSE),KMAX_RESOL=NMAX_RESOL, &
&                 KPROMATR=NPROMATR,KPRGPNS=NPRGPNS,KPRGPEW=NPRGPEW,KPRTRW=NPRTRW, &
&                 KCOMBFLEN=NCOMBFLEN,LDMPOFF=LMPOFF,LDSYNC_TRANS=LSYNC_TRANS, &
&                 LDEQ_REGIONS=LEQ_REGIONS, &
&                 PRAD=ZRA,LDALLOPERM=.TRUE.)

CALL SETUP_TRANS(KSMAX=NSMAX,KDGL=NDGL,KLOEN=NLOEN,LDSPLIT=.TRUE.,&
&                 LDUSEFFTW=.false., LDUSERPNM=LUSERPNM,LDKEEPRPNM=LKEEPRPNM, &
&                 LDUSEFLT=LUSEFLT)
!
CALL TRANS_INQ(KSPEC2=NSPEC2,KSPEC2G=NSPEC2G,KGPTOT=NGPTOT,KGPTOTG=NGPTOTG)

! Default, no blocking
NPROMA=NGPTOT

! Calculate number of NPROMA blocks
NGPBLKS=(NGPTOT-1)/NPROMA+1

!===================================================================================================
! Initialize spectral arrays
!===================================================================================================

! Allocate spectral arrays
! Try to mimick IFS layout as much as possible
NULLIFY(ZVOR)
NULLIFY(ZDIV)
NULLIFY(ZT)
ALLOCATE(SP3D(NFLEVL,NSPEC2,3))
ALLOCATE(ZSP(1,NSPEC2))

SP3D(:,:,:) =  0.0_JPRB
ZSP(:,:)    =  0.0_JPRB

! Initialize all fields to be a randomly chosen spherical harmonic
ZSP(1,162)    = 1.0
SP3D(:,162,:) = 1.0

! Point convenience variables to storage variable SP3D
ZVOR => SP3D(:,:,1)
ZDIV => SP3D(:,:,2)
ZT   => SP3D(:,:,3:3)

!===================================================================================================
! Print information before starting
!===================================================================================================

! PRINT CONFIGURATION DETAILS
WRITE(NOUT,'(A)')'===-=== START OF  RUNTIME PARAMETERS ===-==='
WRITE(NOUT,'(" ")')
WRITE(NOUT,'("NLIN=   ",I10)') NLIN
WRITE(NOUT,'("NQ=     ",I10)') NQ
WRITE(NOUT,'("NSMAX=  ",I10)') NSMAX
WRITE(NOUT,'("NDGL=   ",I10)') NDGL
WRITE(NOUT,'("NPROC=  ",I10)') NPROC
WRITE(NOUT,'("NTHREAD=",I10)') NTHREAD
WRITE(NOUT,'("NPRGPNS=",I10)') NPRGPNS
WRITE(NOUT,'("NPRGPEW=",I10)') NPRGPEW
WRITE(NOUT,'("NPRTRW= ",I10)') NPRTRW
WRITE(NOUT,'("NPRTRV= ",I10)') NPRTRV
WRITE(NOUT,'("NPROMA= ",I10)') NPROMA
WRITE(NOUT,'("NGPTOT= ",I10)') NGPTOT
WRITE(NOUT,'("NGPTOTG=",I10)') NGPTOTG
WRITE(NOUT,'("NFLEVG= ",I10)') NFLEVG
WRITE(NOUT,'("IFLDS=  ",I10)') IFLDS
WRITE(NOUT,'("NSPEC2= ",I10)') NSPEC2
WRITE(NOUT,'("NSPEC2G=",I10)') NSPEC2G
WRITE(NOUT,'("LUSEFLT=",L10)') LUSEFLT
WRITE(NOUT,'(" ")')
WRITE(NOUT,'(A)') '===-=== END OF   RUNTIME PARAMETERS ===-==='

ALLOCATE(IVSET(NFLEVG))

! Compute spectral distribution
ILEV = 0
DO JB=1,NPRTRV
  DO JLEV=1,NUMLL(JB)
    ILEV = ILEV + 1
    IVSET(ILEV) = JB
  ENDDO
ENDDO

ALLOCATE(ITO(IFLDS))
ITO(:)=1

! ALLOCATE GRID-POINT ARRAYS
ALLOCATE(ZWINDS(NPROMA,NFLEVG,4,NGPBLKS))
ALLOCATE(ZGMV(NPROMA,NFLEVG,NDIMGMV,NGPBLKS))
ALLOCATE(ZGMVS(NPROMA,NDIMGMVS,NGPBLKS))

ALLOCATE(ZNORMSP(1))
ALLOCATE(ZNORMSP1(1))
ALLOCATE(ZNORMVOR(NFLEVG))
ALLOCATE(ZNORMVOR1(NFLEVG))
ALLOCATE(ZNORMDIV(NFLEVG))
ALLOCATE(ZNORMDIV1(NFLEVG))
ALLOCATE(ZNORMT(NFLEVG))
ALLOCATE(ZNORMT1(NFLEVG))

IF( VERBOSE ) THEN
  CALL SPECNORM(PSPEC=ZVOR(1:NFLEVL,:),PNORM=ZNORMVOR1,KVSET=IVSET(1:NFLEVG))
  CALL SPECNORM(PSPEC=ZDIV(1:NFLEVL,:),PNORM=ZNORMDIV1,KVSET=IVSET(1:NFLEVG))
  CALL SPECNORM(PSPEC=ZT(1:NFLEVL,:,1),PNORM=ZNORMT1,KVSET=IVSET(1:NFLEVG))
  CALL SPECNORM(PSPEC=ZSP(1:1,:),      PNORM=ZNORMSP1,KVSET=IVSETSC(1:1))

  IF(MYPROC == 1) THEN
    DO IFLD=1,1
      WRITE(NOUT,'("SP  ZNORM(",I4,")=",F20.15)') IFLD,ZNORMSP1(IFLD)
    ENDDO
    DO IFLD=1,NFLEVG
      WRITE(NOUT,'("DIV ZNORM(",I4,")=",F20.15)') IFLD,ZNORMDIV1(IFLD)
    ENDDO
    DO IFLD=1,NFLEVG
      WRITE(NOUT,'("VOR ZNORM(",I4,")=",F20.15)') IFLD,ZNORMVOR1(IFLD)
    ENDDO
    DO IFLD=1,NFLEVG
      WRITE(NOUT,'("T   ZNORM(",I4,")=",F20.15)') IFLD,ZNORMT1(IFLD)
    ENDDO
  ENDIF
ENDIF

ZTINIT=(TIMEF()-ZTINIT)/1000.0_JPRD
WRITE(NOUT,'(" ")')
WRITE(NOUT,'(a,I6,a,F9.2,a)') "TRANSFORM_TEST initialisation, on",NPROC,&
                              & " tasks, took",ZTINIT," sec"
WRITE(NOUT,'(" ")')

IF(ITERS<=0) CALL ABOR1('TRANSFORM_TEST:ITERS <= 0')

ALLOCATE(ZTSTEP(ITERS))
ALLOCATE(ZTSTEP1(ITERS))
ALLOCATE(ZTSTEP2(ITERS))

ZTSTEPAVG=0._JPRD
ZTSTEPMAX=0._JPRD
ZTSTEPMIN=9999999999999999._JPRD
ZTSTEPAVG1=0._JPRD
ZTSTEPMAX1=0._JPRD
ZTSTEPMIN1=9999999999999999._JPRD
ZTSTEPAVG2=0._JPRD
ZTSTEPMAX2=0._JPRD
ZTSTEPMIN2=9999999999999999._JPRD

WRITE(NOUT,'(A)') '===-=== START OF  SPEC TRANSFORMS  ===-==='
WRITE(NOUT,'(" ")')

IF( LSTATS ) THEN
  CALL GSTATS(0,0)
  CALL GSTATS_SETUP(NPROC,MYPROC,NPRCIDS,&
   & LSTATS,LSTATSCPU,LSYNCSTATS,LDETAILED_STATS,LBARRIER_STATS,LBARRIER_STATS2,&
   & LSTATS_OMP,LSTATS_COMMS,LSTATS_MEM,NSTATS_MEM,LSTATS_ALLOC,&
   & LTRACE_STATS,NTRACE_STATS,NPRNT_STATS,LXML_STATS)
  CALL GSTATS_PSUT
  ! TODO: What is this?
  !CALL GSTATS_LABEL_IFS
ENDIF

ZTLOOP=TIMEF()

!===================================================================================================
! Do spectral transform loop
!===================================================================================================

DO JSTEP=1,ITERS
  ZTSTEP(JSTEP)=TIMEF()

  !=================================================================================================
  ! Do inverse transform
  !=================================================================================================

  ZTSTEP1(JSTEP)=TIMEF()
  CALL INV_TRANS(PSPVOR=ZVOR,PSPDIV=ZDIV,PSPSC2=ZSP(1:1,:),&
     & PSPSC3A=ZT,&
     & LDSCDERS=.TRUE.,LDVORGP=.FALSE.,LDDIVGP=.TRUE.,LDUVDER=.FALSE.,&
     & KRESOL=1,KPROMA=NPROMA,KVSETUV=IVSET,KVSETSC2=IVSETSC(1:1),&
     & KVSETSC3A=IVSET,&
     & PGPUV=ZWINDS(:,:,2:4,:),PGP2=ZGMVS(:,1:3,:),&
     & PGP3A=ZGMV(:,:,5:7,:))
  ZTSTEP1(JSTEP)=(TIMEF()-ZTSTEP1(JSTEP))/1000.0_JPRD

  !=================================================================================================
  ! While in grid point space, dump the values to disk
  !=================================================================================================

  ! Dump a field to a binary file
  CALL DUMP_GRIDPOINT_FIELD(JSTEP, MYPROC, NPROMA, NGPBLKS, ZGMVS(:,1,:), 'S', NOUTDUMP)
  CALL DUMP_GRIDPOINT_FIELD(JSTEP, MYPROC, NPROMA, NGPBLKS, ZWINDS(:,NFLEVG,3,:),  'U', NOUTDUMP)
  CALL DUMP_GRIDPOINT_FIELD(JSTEP, MYPROC, NPROMA, NGPBLKS, ZWINDS(:,NFLEVG,4,:),  'V', NOUTDUMP)
  CALL DUMP_GRIDPOINT_FIELD(JSTEP, MYPROC, NPROMA, NGPBLKS, ZGMV(:,NFLEVG,5,:),  'T', NOUTDUMP)

  !=================================================================================================
  ! Do direct transform
  !=================================================================================================

  ZTSTEP2(JSTEP)=TIMEF()
  CALL DIR_TRANS(PSPVOR=ZVOR,PSPDIV=ZDIV,&
      & PSPSC2=ZSP(1:1,:),PSPSC3A=ZT,&
      & KRESOL=1,KPROMA=NPROMA,KVSETUV=IVSET,KVSETSC2=IVSETSC(1:1),&
      & KVSETSC3A=IVSET,&
      & PGPUV=ZWINDS(:,:,3:4,:),PGP2=ZGMVS(:,1:1,:),&
      & PGP3A=ZGMV(:,:,5:5,:))
  ZTSTEP2(JSTEP)=(TIMEF()-ZTSTEP2(JSTEP))/1000.0_JPRD

  !=================================================================================================
  ! Calculate timings
  !=================================================================================================

  ZTSTEP(JSTEP)=(TIMEF()-ZTSTEP(JSTEP))/1000.0_JPRD

  ZTSTEPAVG=ZTSTEPAVG+ZTSTEP(JSTEP)
  ZTSTEPMIN=MIN(ZTSTEP(JSTEP),ZTSTEPMIN)
  ZTSTEPMAX=MAX(ZTSTEP(JSTEP),ZTSTEPMAX)

  ZTSTEPAVG1=ZTSTEPAVG1+ZTSTEP1(JSTEP)
  ZTSTEPMIN1=MIN(ZTSTEP1(JSTEP),ZTSTEPMIN1)
  ZTSTEPMAX1=MAX(ZTSTEP1(JSTEP),ZTSTEPMAX1)

  ZTSTEPAVG2=ZTSTEPAVG2+ZTSTEP2(JSTEP)
  ZTSTEPMIN2=MIN(ZTSTEP2(JSTEP),ZTSTEPMIN2)
  ZTSTEPMAX2=MAX(ZTSTEP2(JSTEP),ZTSTEPMAX2)

  !=================================================================================================
  ! Print norms
  !=================================================================================================

  IF( VERBOSE )THEN
    CALL SPECNORM(PSPEC=ZSP(1:1,:),       PNORM=ZNORMSP, KVSET=IVSETSC(1:1))
    CALL SPECNORM(PSPEC=ZVOR(1:NFLEVL,:), PNORM=ZNORMVOR,KVSET=IVSET(1:NFLEVG))
    CALL SPECNORM(PSPEC=ZDIV(1:NFLEVL,:), PNORM=ZNORMDIV,KVSET=IVSET(1:NFLEVG))
    CALL SPECNORM(PSPEC=ZT(1:NFLEVL,:,1), PNORM=ZNORMT, KVSET=IVSET(1:NFLEVG))

    IF( MYPROC==1 ) THEN
      ! SURFACE PRESSURE
      ZMAXERR(:)=-999.0
      DO IFLD=1,1
        ZERR(1)=ABS(ZNORMSP1(IFLD)/ZNORMSP(IFLD)-1.0_JPRB)
        ZMAXERR(1)=MAX(ZMAXERR(1),ZERR(1))
      ENDDO
      ! DIVERGENCE
      DO IFLD=1,NFLEVG
        ZERR(2)=ABS(ZNORMDIV1(IFLD)/ZNORMDIV(IFLD)-1.0_JPRB)
        ZMAXERR(2)=MAX(ZMAXERR(2),ZERR(2))
      ENDDO
      ! VORTICITY
      DO IFLD=1,NFLEVG
        ZERR(3)=ABS(ZNORMVOR1(IFLD)/ZNORMVOR(IFLD)-1.0_JPRB)
        ZMAXERR(3)=MAX(ZMAXERR(3),ZERR(3))
      ENDDO
      ! TEMPERATURE
      DO IFLD=1,NFLEVG
        ZERR(4)=ABS(ZNORMT1(IFLD)/ZNORMT(IFLD)-1.0_JPRB)
        ZMAXERR(4)=MAX(ZMAXERR(4),ZERR(4))
      ENDDO
      WRITE(NOUT,'("time step ",I6," took", F8.4," | SP max err="E10.3,&
                 & " | DIV max err="E10.3," | VOR max err="E10.3," | T max err="E10.3)') &
                 &  JSTEP,ZTSTEP(JSTEP),ZMAXERR(1),ZMAXERR(2),ZMAXERR(3),ZMAXERR(4)
    ENDIF
  ELSE
    WRITE(NOUT,'("time step ",I6," took", F8.4)') JSTEP,ZTSTEP(JSTEP)
  ENDIF
ENDDO

ZTLOOP=(TIMEF()-ZTLOOP)/1000.0_JPRD

WRITE(NOUT,'(" ")')
WRITE(NOUT,'(A)') '===-=== END OF   SPEC TRANSFORMS  ===-==='
WRITE(NOUT,'(" ")')


IF( VERBOSE ) THEN
  CALL SPECNORM(PSPEC=ZVOR(1:NFLEVL,:),PNORM=ZNORMVOR,KVSET=IVSET(1:NFLEVG))
  CALL SPECNORM(PSPEC=ZDIV(1:NFLEVL,:),PNORM=ZNORMDIV,KVSET=IVSET(1:NFLEVG))
  CALL SPECNORM(PSPEC=ZT(1:NFLEVL,:,1),PNORM=ZNORMT,KVSET=IVSET(1:NFLEVG))
  CALL SPECNORM(PSPEC=ZSP(1:1,:),      PNORM=ZNORMSP,KVSET=IVSETSC(1:1))

  IF(MYPROC == 1) THEN
    ! SURFACE PRESSURE
    ZMAXERR(:)=-999.0
    DO IFLD=1,1
      ZERR(1)=ABS(ZNORMSP1(IFLD)/ZNORMSP(IFLD)-1.0D0)
      ZMAXERR(1)=MAX(ZMAXERR(1),ZERR(1))
      WRITE(NOUT,'("SP ZNORM(",I4,")=",F20.15," ERR=",E10.3)') IFLD,ZNORMSP(IFLD),ZERR(1)
    ENDDO
    ! DIVERGENCE
    DO IFLD=1,NFLEVG
      ZERR(2)=ABS(ZNORMDIV1(IFLD)/ZNORMDIV(IFLD)-1.0D0)
      ZMAXERR(2)=MAX(ZMAXERR(2),ZERR(2))
      WRITE(NOUT,'("DIV ZNORM(",I4,")=",F20.15," ERR=",E10.3)') IFLD,ZNORMDIV(IFLD),ZERR(2)
    ENDDO
    ! VORTICITY
    DO IFLD=1,NFLEVG
      ZERR(3)=ABS(ZNORMVOR1(IFLD)/ZNORMVOR(IFLD)-1.0D0)
      ZMAXERR(3)=MAX(ZMAXERR(3),ZERR(3))
      WRITE(NOUT,'("VOR ZNORM(",I4,")=",F20.15," ERR=",E10.3)') IFLD,ZNORMVOR(IFLD),ZERR(3)
    ENDDO
    ! TEMPERATURE
    DO IFLD=1,NFLEVG
      ZERR(4)=ABS(ZNORMT1(IFLD)/ZNORMT(IFLD)-1.0D0)
      ZMAXERR(4)=MAX(ZMAXERR(4),ZERR(4))
      WRITE(NOUT,'("T ZNORM(",I4,")=",F20.15," ERR=",E10.3)') IFLD,ZNORMT(IFLD),ZERR(4)
    ENDDO
    ! MAXIMUM ERROR ACROSS ALL FIELDS
    ZMAXERRG=MAX(MAX(ZMAXERR(1),ZMAXERR(2)),MAX(ZMAXERR(2),ZMAXERR(3)))

    WRITE(NOUT,'("SURFACE PRESSURE MAX ERROR=",E10.3)')ZMAXERR(1)
    WRITE(NOUT,'("DIVERGENCE       MAX ERROR=",E10.3)')ZMAXERR(2)
    WRITE(NOUT,'("VORTICITY        MAX ERROR=",E10.3)')ZMAXERR(3)
    WRITE(NOUT,'("TEMPERATURE      MAX ERROR=",E10.3)')ZMAXERR(4)
    WRITE(NOUT,'("GLOBAL           MAX ERROR=",E10.3)')ZMAXERRG

  ENDIF
ENDIF

CALL MPL_ALLREDUCE(ZTLOOP,     'SUM', LDREPROD=.FALSE.)
CALL MPL_ALLREDUCE(ZTSTEP,     'SUM', LDREPROD=.FALSE.)
CALL MPL_ALLREDUCE(ZTSTEPAVG,  'SUM', LDREPROD=.FALSE.)
CALL MPL_ALLREDUCE(ZTSTEPMAX,  'MAX', LDREPROD=.FALSE.)
CALL MPL_ALLREDUCE(ZTSTEPMIN,  'MIN', LDREPROD=.FALSE.)

CALL MPL_ALLREDUCE(ZTSTEP1,    'SUM', LDREPROD=.FALSE.)
CALL MPL_ALLREDUCE(ZTSTEPAVG1, 'SUM', LDREPROD=.FALSE.)
CALL MPL_ALLREDUCE(ZTSTEPMAX1, 'MAX', LDREPROD=.FALSE.)
CALL MPL_ALLREDUCE(ZTSTEPMIN1, 'MIN', LDREPROD=.FALSE.)

CALL MPL_ALLREDUCE(ZTSTEP2,    'SUM', LDREPROD=.FALSE.)
CALL MPL_ALLREDUCE(ZTSTEPAVG2, 'SUM', LDREPROD=.FALSE.)
CALL MPL_ALLREDUCE(ZTSTEPMAX2, 'MAX', LDREPROD=.FALSE.)
CALL MPL_ALLREDUCE(ZTSTEPMIN2, 'MIN', LDREPROD=.FALSE.)


ZTSTEPAVG=(ZTSTEPAVG/REAL(NPROC,JPRB))/REAL(ITERS,JPRD)
ZTLOOP=ZTLOOP/REAL(NPROC,JPRD)
ZTSTEP(:)=ZTSTEP(:)/REAL(NPROC,JPRD)

CALL SORT(ZTSTEP,ITERS)
ZTSTEPMED = ZTSTEP(ITERS/2)

ZTSTEPAVG1=(ZTSTEPAVG1/REAL(NPROC,JPRB))/REAL(ITERS,JPRD)
ZTSTEP1(:)=ZTSTEP1(:)/REAL(NPROC,JPRD)

CALL SORT(ZTSTEP1,ITERS)
ZTSTEPMED1 = ZTSTEP1(ITERS/2)

ZTSTEPAVG2=(ZTSTEPAVG2/REAL(NPROC,JPRB))/REAL(ITERS,JPRD)
ZTSTEP2(:)=ZTSTEP2(:)/REAL(NPROC,JPRD)

CALL SORT(ZTSTEP2,ITERS)
ZTSTEPMED2 = ZTSTEP2(ITERS/2)

IF(MYPROC == 1)THEN
  WRITE(NOUT,'(" ")')
  WRITE(NOUT,'(A)') '===-=== START   OF TIME STEP STATS ===-==='
  WRITE(NOUT,'(" ")')
  WRITE(NOUT,'("INVERSE TRANSFORMS")')
  WRITE(NOUT,'("------------------")')
  WRITE(NOUT,'("AVG  (s): ",F8.4)') ZTSTEPAVG1
  WRITE(NOUT,'("MIN  (s): ",F8.4)') ZTSTEPMIN1
  WRITE(NOUT,'("MAX  (s): ",F8.4)') ZTSTEPMAX1
  WRITE(NOUT,'("MED  (s): ",F8.4)') ZTSTEPMED1
  WRITE(NOUT,'(" ")')
  WRITE(NOUT,'("DIRECT TRANSFORMS")')
  WRITE(NOUT,'("-----------------")')
  WRITE(NOUT,'("AVG  (s): ",F8.4)') ZTSTEPAVG2
  WRITE(NOUT,'("MIN  (s): ",F8.4)') ZTSTEPMIN2
  WRITE(NOUT,'("MAX  (s): ",F8.4)') ZTSTEPMAX2
  WRITE(NOUT,'("MED  (s): ",F8.4)') ZTSTEPMED2
  WRITE(NOUT,'(" ")')
  WRITE(NOUT,'("INVERSE-DIRECT TRANSFORMS")')
  WRITE(NOUT,'("-------------------------")')
  WRITE(NOUT,'("AVG  (s): ",F8.4)') ZTSTEPAVG
  WRITE(NOUT,'("MIN  (s): ",F8.4)') ZTSTEPMIN
  WRITE(NOUT,'("MAX  (s): ",F8.4)') ZTSTEPMAX
  WRITE(NOUT,'("MED  (s): ",F8.4)') ZTSTEPMED
  WRITE(NOUT,'("LOOP (s): ",F8.4)') ZTLOOP
  WRITE(NOUT,'(" ")')
  WRITE(NOUT,'(A)') '===-=== END     OF TIME STEP STATS ===-==='
  WRITE(NOUT,'(" ")')
ENDIF

IF( LSTACK ) THEN
!           gather stack usage statistics
  ISTACK = GETSTACKUSAGE()
  IF(MYPROC == 1) THEN
    PRINT 9000, istack
    9000 FORMAT("Stack Utilisation Information",/,&
         &"=============================",//,&
         &"Task           Size(Bytes)",/,&
         &"====           ===========",//,&
         &"   1",11x,I10)

    DO I=2,NPROC
      CALL MPL_RECV(ISTACK,KSOURCE=NPRCIDS(I),KTAG=I, &
           & CDSTRING='TRANSFORM_TEST:')
      PRINT '(I4,11X,I10)', I,ISTACK
    ENDDO
  ELSE
    CALL MPL_SEND(ISTACK,KDEST=NPRCIDS(1),KTAG=MYPROC, &
             &   CDSTRING='TRANSFORM_TEST:')
  ENDIF
ENDIF


!===================================================================================================

IF( LSTATS ) THEN
  CALL GSTATS(0,1)
  CALL GSTATS_PRINT(NOUT,ZAVEAVE,JPMAXSTAT)
ENDIF

!===================================================================================================

! CLOSE FILE
IF( NPROC > 1 ) THEN
  IF( MYPROC /= 1 ) THEN
    CLOSE(UNIT=NOUT)
  ENDIF
ENDIF

DEALLOCATE(ZWINDS)
DEALLOCATE(ZGMV)
DEALLOCATE(ZGMVS)

!===================================================================================================
! Finalize MPI
!===================================================================================================

CALL MPL_BARRIER()
CALL MPL_END()

!===================================================================================================

CONTAINS

!===================================================================================================

SUBROUTINE SORT(A, N)
    IMPLICIT NONE
    INTEGER(KIND=JPIM) :: N, I, J
    REAL(KIND=JPRD)    :: A(N), X

    DO I = 2, N
        X = A(I)
        J = I - 1
        DO WHILE (J >= 1)
            IF (A(J) <= X) EXIT
            A(J + 1) = A(J)
            J = J - 1
        END DO
        A(J + 1) = X
    END DO
END SUBROUTINE

!===================================================================================================

SUBROUTINE DUMP_GRIDPOINT_FIELD(JSTEP, MYPROC, NPROMA, NGPBLKS, FLD, FLDCHAR, NOUTDUMP)

! Dump a 2D field to a binary file.

INTEGER(KIND=JPIM), INTENT(IN) :: JSTEP ! Time step, used for naming file
INTEGER(KIND=JPIM), INTENT(IN) :: MYPROC ! MPI rank, used for naming file
INTEGER(KIND=JPIM), INTENT(IN) :: NPROMA ! Size of NPROMA
INTEGER(KIND=JPIM), INTENT(IN) :: NGPBLKS ! Number of NPROMA blocks
REAL(KIND=JPRB)   , INTENT(IN) :: FLD(NPROMA,NGPBLKS) ! 2D field
CHARACTER         , INTENT(IN) :: FLDCHAR ! Single character field identifier
INTEGER(KIND=JPIM), INTENT(IN) :: NOUTDUMP ! Unit number for output file

CHARACTER(LEN=14) :: FILENAME = "X.XXX.XXXX.dat"

WRITE(FILENAME(1:1),'(A1)') FLDCHAR
WRITE(FILENAME(3:5),'(I3.3)') JSTEP
WRITE(FILENAME(7:10),'(I4.4)') MYPROC

OPEN(NOUTDUMP, FILE=FILENAME, FORM="UNFORMATTED")
WRITE(NOUTDUMP) RESHAPE(FLD, (/ NPROMA*NGPBLKS /))
CLOSE(NOUTDUMP)

END SUBROUTINE DUMP_GRIDPOINT_FIELD

END PROGRAM TRANSFORM_TEST

!===================================================================================================
