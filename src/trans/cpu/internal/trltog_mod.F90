! (C) Copyright 1995- ECMWF.
! (C) Copyright 1995- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TRLTOG_MOD

USE PARKIND1, ONLY : JPIM

IMPLICIT NONE

PRIVATE TRLTOG_PROLOG, TRLTOG_COMM

CONTAINS

SUBROUTINE TRLTOG(PGLAT,KF_FS,KF_GP,KF_SCALARS_G,KVSET,KPTRGP,&
 &PGP,PGPUV,PGP3A,PGP3B,PGP2)

!**** *TRLTOG * - head routine for transposition of grid point data from latitudinal
!                 to column structure (this takes place between inverse
!                 FFT and grid point calculations)
!                 TRLTOG is the inverse of TRGTOL

!**   Interface.
!     ----------
!        *call* *TRLTOG(...)

!        Explicit arguments :
!        --------------------
!           PGLAT    -  Latitudinal data ready for direct FFT (input)
!           PGP    -  Blocked grid point data    (output)
!           KVSET    - "v-set" for each field      (input)

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        R. El Khatib *Meteo-France*

!     Modifications.
!     --------------
!        Original  : 18-Aug-2014 from trltog
!        R. El Khatib 09-Sep-2020 NSTACK_MEMORY_TR
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DISTR       ,ONLY : D
USE TRGL_MOD, ONLY : TRGL_CTX, TRGL_VARS
IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KF_FS,KF_GP
INTEGER(KIND=JPIM),INTENT(IN)  :: KF_SCALARS_G
REAL(KIND=JPRB),INTENT(IN)     :: PGLAT(KF_FS,D%NLENGTF)
INTEGER(KIND=JPIM), INTENT(IN) :: KVSET(KF_GP)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KPTRGP(:)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT)     :: PGP(:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT)     :: PGPUV(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT)     :: PGP3A(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT)     :: PGP3B(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT)     :: PGP2(:,:,:)

TYPE(TRGL_CTX) :: YLCTX

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TRLTOG',0,ZHOOK_HANDLE)

CALL TRLTOG_PROLOG(KF_FS,KF_GP,KVSET,YLCTX)

CALL TRLTOG_COMM(PGLAT,KF_FS,KF_GP,KF_SCALARS_G,KVSET, &
                 & KPTRGP,PGP,PGPUV,PGP3A,PGP3B,PGP2,YLCTX)

IF (LHOOK) CALL DR_HOOK('TRLTOG',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE TRLTOG

SUBROUTINE TRLTOG_PROLOG(KF_FS,KF_GP,KVSET,YDCTX)
 
!**** *TRLTOG_PROLOG * - prolog for transposition of grid point data from latitudinal
!                 to column structure (this takes place between inverse
!                 FFT and grid point calculations) : the purpose is essentially 
!                 to compute the size of communication buffers in order to enable
!                 the use of automatic arrays later.
!                 TRLTOG_PROLOG is the inverse of TRGTOL_PROLOG

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *call* *TRLTOG_PROLOG(...)

!        Explicit arguments :
!        --------------------
!           KVSET    - "v-set" for each field      (input)

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        R. El Khatib *Meteo-France*

!     Modifications.
!     --------------
!        Original  : 18-Aug-2014 from trltog
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM

USE INIGPTR_MOD     ,ONLY : INIGPTR
USE PE2SET_MOD      ,ONLY : PE2SET
USE TRGL_MOD  ,ONLY : TRGL_CTX, ALLOCATE_CTX, TRGL_PROLOG

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KF_FS,KF_GP
TYPE (TRGL_CTX), INTENT(INOUT) :: YDCTX
INTEGER(KIND=JPIM),INTENT(IN) :: KVSET(KF_GP)

!     ------------------------------------------------------------------
!*       0.    Some initializations
!              --------------------

CALL ALLOCATE_CTX(YDCTX, KF_GP)
CALL GSTATS(1806,0)
CALL TRGL_PROLOG(KF_FS,KF_GP,KVSET,YDCTX, .TRUE.)  ! TRUE. for LSEND
CALL GSTATS(1806,1)

END SUBROUTINE TRLTOG_PROLOG


SUBROUTINE TRLTOG_COMM(PGLAT,KF_FS,KF_GP,KF_SCALARS_G,KVSET,&
 & KPTRGP,PGP,PGPUV,PGP3A,PGP3B,PGP2,YDCTX)
 

!**** *trltog * - transposition of grid point data from latitudinal
!                 to column structure. This takes place between inverse
!                 FFT and grid point calculations.
!                 TRLTOG_COMM is the inverse of TRGTOL

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *call* *trltog(...)

!        Explicit arguments :
!        --------------------
!           PGLAT    -  Latitudinal data ready for direct FFT (input)
!           PGP    -  Blocked grid point data    (output)
!           KVSET    - "v-set" for each field      (input)

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        MPP Group *ECMWF*

!     Modifications.
!     --------------
!        Original  : 95-10-01
!        D.Dent    : 97-08-04 Reorganisation to allow NPRTRV
!                             to differ from NPRGPEW
!        =99-03-29= Mats Hamrud and Deborah Salmond
!                   JUMP in FFT's changed to 1
!                   KINDEX introduced and PCOMBUF not used for same PE
!         01-11-23  Deborah Salmond and John Hague
!                   LIMP_NOOLAP Option for non-overlapping message passing
!                               and buffer packing
!         01-12-18  Peter Towers
!                   Improved vector performance of LTOG_PACK,LTOG_UNPACK
!         03-0-02   G. Radnoti: Call barrier always when nproc>1
!         08-01-01  G.Mozdzynski: cleanup
!         09-01-02  G.Mozdzynski: use non-blocking recv and send
!        R. El Khatib 09-Sep-2020 64 bits addressing for PGLAT
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB    ,JPIB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE MPL_MODULE  ,ONLY : MPL_RECV, MPL_SEND, MPL_WAIT, JP_NON_BLOCKING_STANDARD, MPL_WAITANY, &
     & JP_BLOCKING_STANDARD, MPL_BARRIER, JP_BLOCKING_BUFFERED

USE TPM_GEN         ,ONLY : NTRANS_SYNC_LEVEL, NSTACK_MEMORY_TR
USE TPM_DISTR       ,ONLY : D, MTAGLG,      &
     &                      NPRCIDS, MYPROC, NPROC

USE PE2SET_MOD      ,ONLY : PE2SET
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
USE TRGL_MOD  ,ONLY : TRGL_CTX, TRGL_VARS, &
                             &TRGL_ALLOCATE_VARS, TRGL_ALLOCATE_HEAP_BUFFER,&
                             &TRGL_INIT_VARS, TRGL_INIT_OFF_VARS, &
                             &TGRL_COPY_ZCOMBUF, TGRL_COPY_PGLAT, TGRL_INIT_PACKING_VARS

IMPLICIT NONE


INTEGER(KIND=JPIM), INTENT(IN) :: KF_FS,KF_GP
REAL(KIND=JPRB),INTENT(IN)     :: PGLAT(KF_FS,D%NLENGTF)
INTEGER(KIND=JPIM), INTENT(IN) :: KVSET(KF_GP)
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCALARS_G
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KPTRGP(:)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT) :: PGP(:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT) :: PGPUV(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT) :: PGP3A(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT) :: PGP3B(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT) :: PGP2(:,:,:)

TYPE (TRGL_CTX), INTENT(INOUT), TARGET :: YDCTX
! LOCAL VARIABLES
TYPE(TRGL_VARS) :: YLVARS

INTEGER(KIND=JPIM) :: IREQ_SEND(NPROC)
INTEGER(KIND=JPIM) :: IREQ_RECV(NPROC)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: IRECV
INTEGER(KIND=JPIM) :: ISEND, ITAG, JL, JFLD, INS, INR, JNR
INTEGER(KIND=JPIM) :: II,ILEN
INTEGER(KIND=JPIM) :: ISEND_FLD_START,ISEND_FLD_END

!     LOCAL ARRAYS
REAL(KIND=JPRB), TARGET :: ZCOMBUFS_STACK(-1:YDCTX%ISENDCOUNT,MERGE (YDCTX%INSEND,0,NSTACK_MEMORY_TR/=0))
REAL(KIND=JPRB), TARGET :: ZCOMBUFR_STACK(-1:YDCTX%IRECVCOUNT,MERGE (YDCTX%INRECV,0,NSTACK_MEMORY_TR/=0))

REAL(KIND=JPRB), ALLOCATABLE, TARGET, SAVE :: ZCOMBUFS_HEAP(:,:)
REAL(KIND=JPRB), ALLOCATABLE, TARGET, SAVE :: ZCOMBUFR_HEAP(:,:)

REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZCOMBUFS(:,:)
REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZCOMBUFR(:,:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE_BAR

!     ------------------------------------------------------------------

!*       0.    Some initializations
!              --------------------
ASSOCIATE (KSENDCOUNT=>YDCTX%ISENDCOUNT, KRECVCOUNT=>YDCTX%IRECVCOUNT, KNSEND=>YDCTX%INSEND, KNRECV=>YDCTX%INRECV, &
         & KSENDTOT=>YDCTX%ISENDTOT, KRECVTOT=>YDCTX%IRECVTOT, KSEND=>YDCTX%ISEND, KRECV=>YDCTX%IRECV, KINDEX=>YDCTX%IINDEX, &
         & KNDOFF=>YDCTX%INDOFF, KGPTRSEND =>YDCTX%IGPTRSEND, KSETAL=>YDCTX%ISETAL, KSETBL=>YDCTX%ISETBL, KSETWL=>YDCTX%ISETWL)

IF (NSTACK_MEMORY_TR == 0) THEN
  CALL TRGL_ALLOCATE_HEAP_BUFFER(ZCOMBUFS_HEAP, YDCTX%ISENDCOUNT, YDCTX%INSEND)
  CALL TRGL_ALLOCATE_HEAP_BUFFER(ZCOMBUFR_HEAP, YDCTX%IRECVCOUNT, YDCTX%INRECV)

! Now, force the OS to allocate this shared array right now, not when it starts to be used which is
! an OPEN-MP loop, that would cause a threads synchronization lock :
  IF (YDCTX%INSEND > 0 .AND. YDCTX%ISENDCOUNT >=-1) ZCOMBUFS_HEAP(-1,1)=HUGE(1._JPRB)
  ZCOMBUFS (-1:,1:) => ZCOMBUFS_HEAP
  ZCOMBUFR (-1:,1:) => ZCOMBUFR_HEAP
ELSE
  ZCOMBUFS (-1:,1:) => ZCOMBUFS_STACK
  ZCOMBUFR (-1:,1:) => ZCOMBUFR_STACK
ENDIF

ITAG = MTAGLG

IF (LHOOK) CALL DR_HOOK('TRLTOG_BAR',0,ZHOOK_HANDLE_BAR)
CALL GSTATS_BARRIER(762)
IF (LHOOK) CALL DR_HOOK('TRLTOG_BAR',1,ZHOOK_HANDLE_BAR)

CALL GSTATS(805,0)

IF (NTRANS_SYNC_LEVEL <= 0) THEN
  !...Receive loop.........................................................
  DO INR=1,KNRECV
    IRECV=KRECV(INR)
    CALL MPL_RECV(ZCOMBUFR(-1:KRECVTOT(IRECV),INR), &
           & KSOURCE=NPRCIDS(IRECV), &
           & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQ_RECV(INR), &
           & KTAG=ITAG,CDSTRING='TRLTOG_COMM: NON-BLOCKING IRECV' )
  ENDDO
ENDIF

CALL GSTATS(805,1)

CALL GSTATS(1806,0)
YDCTX%LLINDER = .FALSE.
YDCTX%LLPGPONLY = .FALSE.
IF(PRESENT(KPTRGP))  YDCTX%LLINDER = .TRUE.
IF(PRESENT(PGP))     YDCTX%LLPGPONLY = .TRUE.
CALL TRGL_ALLOCATE_VARS(YLVARS, KF_GP)
CALL TRGL_INIT_VARS(YLVARS, KF_SCALARS_G, PGP, PGPUV, PGP3A, PGP3B, PGP2)

CALL GSTATS(1806,1)


! Copy local contribution

IF(KSENDTOT(MYPROC) > 0 )THEN
  CALL TRGL_INIT_OFF_VARS(YDCTX,YLVARS,KVSET,KPTRGP,KF_GP)
  CALL GSTATS(1604,0)
  CALL TGRL_COPY_PGLAT(PGLAT, YDCTX, YLVARS, PGP, PGPUV,PGP3A, PGP3B,PGP2)
  CALL GSTATS(1604,1)
ENDIF

!
! loop over the number of processors we need to communicate with.
! NOT MYPROC
!
! Now overlapping buffer packing/unpacking with sends/waits
! Time as if all communications to avoid double accounting

CALL GSTATS(805,0)

!  Pack+send loop.........................................................

ISEND_FLD_START = 1
ISEND_FLD_END   = KF_FS
DO INS=1,KNSEND
  ISEND=KSEND(INS)
  ILEN = KSENDTOT(ISEND)/KF_FS
  !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JFLD,JL,II)
  DO JL=1,ILEN
    II = KINDEX(KNDOFF(ISEND)+JL)
    DO JFLD=ISEND_FLD_START,ISEND_FLD_END
      ZCOMBUFS((JFLD-ISEND_FLD_START)*ILEN+JL,INS) = PGLAT(JFLD,II)
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  ZCOMBUFS(-1,INS) = 1
  ZCOMBUFS(0,INS)  = KF_FS
  IF (NTRANS_SYNC_LEVEL <= 1) THEN
    CALL MPL_SEND(ZCOMBUFS(-1:KSENDTOT(ISEND),INS),KDEST=NPRCIDS(ISEND),&
          & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQ_SEND(INS), &
          & KTAG=ITAG,CDSTRING='TRLTOG_COMM: NON-BLOCKING ISEND')
  ELSE
    CALL MPL_SEND(ZCOMBUFS(-1:KSENDTOT(ISEND),INS),KDEST=NPRCIDS(ISEND),&
          & KMP_TYPE=JP_BLOCKING_BUFFERED, &
          & KTAG=ITAG,CDSTRING='TRLTOG_COMM: BLOCKING BUFFERED BSEND')
  ENDIF
ENDDO

!  Unpack loop.........................................................

CALL TGRL_INIT_PACKING_VARS(.FALSE., YDCTX,YLVARS, KVSET, KF_GP)

DO JNR=1,KNRECV

  IF (NTRANS_SYNC_LEVEL <= 0) THEN
    CALL MPL_WAITANY(KREQUEST=IREQ_RECV(1:KNRECV),KINDEX=INR,&
          & CDSTRING='TRLTOG_COMM: WAIT FOR ANY RECEIVES')
  ELSE
    INR = JNR
    IRECV=KRECV(INR)
    CALL MPL_RECV(ZCOMBUFR(-1:KRECVTOT(IRECV),INR), &
          & KSOURCE=NPRCIDS(IRECV), &
          & KMP_TYPE=JP_BLOCKING_STANDARD, &
          & KTAG=ITAG,CDSTRING='TRLTOG_COMM: BLOCKING RECV' )
  ENDIF

  CALL TGRL_COPY_ZCOMBUF(YDCTX, YLVARS, INR, &
                      &  ZCOMBUFR, &
                      & KPTRGP, &
                      & PGP, PGPUV, PGP3A, PGP3B, PGP2) 
ENDDO

IF (NTRANS_SYNC_LEVEL <= 1) THEN
  IF(KNSEND > 0) THEN
    CALL MPL_WAIT(KREQUEST=IREQ_SEND(1:KNSEND),CDSTRING='TRLTOG_COMM: WAIT FOR ISENDS')
  ENDIF
ENDIF

IF (NTRANS_SYNC_LEVEL >= 1) THEN
  CALL MPL_BARRIER(CDSTRING='TRLTOG_COMM: BARRIER AT END')
ENDIF

CALL GSTATS(805,1)

CALL GSTATS_BARRIER2(762)

END ASSOCIATE

END SUBROUTINE TRLTOG_COMM
END MODULE TRLTOG_MOD
