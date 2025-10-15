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
USE TRXTOY_CTX_MOD, ONLY : TRXTOY_CTX, TRXTOY_VARS
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

TYPE(TRXTOY_CTX) :: YLCTX

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

USE TPM_DISTR       ,ONLY : D, MYSETW, NPRTRNS, MYPROC, NPROC
USE TPM_TRANS       ,ONLY : NGPBLKS

USE INIGPTR_MOD     ,ONLY : INIGPTR
USE PE2SET_MOD      ,ONLY : PE2SET
USE TRXTOY_CTX_MOD  ,ONLY : TRXTOY_CTX, ALLOCATE_CTX

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KF_FS,KF_GP
TYPE (TRXTOY_CTX), INTENT(INOUT) :: YDCTX
INTEGER(KIND=JPIM),INTENT(IN) :: KVSET(KF_GP)

INTEGER(KIND=JPIM) :: IGPTRRECV(NPRTRNS)
INTEGER(KIND=JPIM) :: IFIRSTLAT, IGL, IGLL, ILASTLAT, IPOS, ISETA, ISETB, ISETV
INTEGER(KIND=JPIM) :: JFLD, JGL, JL, ISETW, JROC, J
INTEGER(KIND=JPIM) :: INDOFFX

!     ------------------------------------------------------------------
!*       0.    Some initializations
!              --------------------


CALL ALLOCATE_CTX(YDCTX, NPROC,D%NLENGTF,NGPBLKS,NPRTRNS)

CALL GSTATS(1806,0)

CALL INIGPTR(YDCTX%IGPTRSEND,IGPTRRECV)

INDOFFX  = 0

YDCTX%INRECV = 0
YDCTX%INSEND = 0

YDCTX%INDOFF = -999

DO JROC=1,NPROC

  CALL PE2SET(JROC,YDCTX%ISETAL(JROC),YDCTX%ISETBL(JROC),YDCTX%ISETWL(JROC),YDCTX%ISETVL(JROC))
  
  ISETA=YDCTX%ISETAL(JROC)
  ISETB=YDCTX%ISETBL(JROC)
  ISETW=YDCTX%ISETWL(JROC)
  ISETV=YDCTX%ISETVL(JROC)
!             count up expected number of fields
  IPOS = 0
  DO JFLD=1,KF_GP
    IF(KVSET(JFLD) == ISETV .OR. KVSET(JFLD) == -1) IPOS = IPOS+1
  ENDDO
  YDCTX%IRECVTOT(JROC) = IGPTRRECV(ISETW)*IPOS
  IF(YDCTX%IRECVTOT(JROC) > 0 .AND. MYPROC /= JROC) THEN
    YDCTX%INRECV = YDCTX%INRECV + 1
    YDCTX%IRECV(YDCTX%INRECV)=JROC
  ENDIF

  IFIRSTLAT = MAX(D%NPTRLS(MYSETW),D%NFRSTLAT(ISETA))
  ILASTLAT  = MIN(D%NPTRLS(MYSETW)+D%NULTPP(MYSETW)-1,D%NLSTLAT(ISETA))

  IPOS = 0
  DO JGL=IFIRSTLAT,ILASTLAT
    IGL  = D%NPTRFRSTLAT(ISETA)+JGL-D%NFRSTLAT(ISETA)
    IPOS = IPOS+D%NONL(IGL,ISETB)
  ENDDO

  YDCTX%ISENDTOT(JROC) = IPOS*KF_FS
  IF( JROC /= MYPROC) THEN    
    IF(YDCTX%ISENDTOT(JROC) > 0) THEN
      YDCTX%INSEND = YDCTX%INSEND+1
      YDCTX%ISEND(YDCTX%INSEND)=JROC
    ENDIF
  ENDIF

  IF(IPOS > 0) THEN
    YDCTX%INDOFF(JROC) = INDOFFX
    INDOFFX = INDOFFX+IPOS
    IPOS = 0
    DO JGL=IFIRSTLAT,ILASTLAT
      IGL  = D%NPTRFRSTLAT(ISETA)+JGL-D%NFRSTLAT(ISETA)
      IGLL = JGL-D%NPTRLS(MYSETW)+1
      DO JL=D%NSTA(IGL,ISETB)+D%NSTAGTF(IGLL),&
       &D%NSTA(IGL,ISETB)+D%NSTAGTF(IGLL)+D%NONL(IGL,ISETB)-1
        IPOS = IPOS+1
        YDCTX%IINDEX(IPOS+YDCTX%INDOFF(JROC)) = JL
      ENDDO
    ENDDO
  ENDIF
ENDDO

YDCTX%ISENDCOUNT=0
YDCTX%IRECVCOUNT=0
DO J=1,NPROC
  YDCTX%ISENDCOUNT=MAX(YDCTX%ISENDCOUNT,YDCTX%ISENDTOT(J))
  YDCTX%IRECVCOUNT=MAX(YDCTX%IRECVCOUNT,YDCTX%IRECVTOT(J))
ENDDO

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
USE TPM_TRANS       ,ONLY : NGPBLKS

USE PE2SET_MOD      ,ONLY : PE2SET
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
USE TRXTOY_CTX_MOD  ,ONLY : TRXTOY_CTX, TRXTOY_VARS, &
                             &ALLOCATE_TOG_VARS, ALLOCATE_HEAP_BUFFER,&
                             &INIT_TOG_VARS, INIT_TOG_OFF_VARS, &
                             &COPY_ZCOMBUF, COPY_PGLAT

IMPLICIT NONE


INTEGER(KIND=JPIM), INTENT(IN) :: KF_FS,KF_GP
REAL(KIND=JPRB),INTENT(IN)     :: PGLAT(KF_FS,D%NLENGTF)
INTEGER(KIND=JPIM), INTENT(IN) :: KVSET(KF_GP)
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCALARS_G
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KPTRGP(:)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT)     :: PGP(:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT)     :: PGPUV(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT)     :: PGP3A(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT)     :: PGP3B(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT)     :: PGP2(:,:,:)

TYPE (TRXTOY_CTX), INTENT(INOUT), TARGET :: YDCTX
! LOCAL VARIABLES
TYPE(TRXTOY_VARS) :: YDTOGLVARS
INTEGER(KIND=JPIM) :: IPOSPLUS(YDCTX%INRECV)
INTEGER(KIND=JPIM) :: ISETW(YDCTX%INRECV)
INTEGER(KIND=JPIM) :: IJPOS(NGPBLKS,YDCTX%INRECV)
INTEGER(KIND=JPIM) :: IFLDA(KF_GP,YDCTX%INRECV)
INTEGER(KIND=JPIM) :: IREQ_SEND(NPROC)
INTEGER(KIND=JPIM) :: IREQ_RECV(NPROC)

INTEGER(KIND=JPIM) :: IFIRST, IFLD, ILAST, IPOS, IRECV, ISETV
INTEGER(KIND=JPIM) :: ISEND, ITAG,  JBLK, JFLD, JL, IFLDS, INR, INS
INTEGER(KIND=JPIM) :: II,ILEN

LOGICAL   :: LLPGPONLY, LLINDER
INTEGER(KIND=JPIM) :: JNR
INTEGER(KIND=JPIM) :: ISEND_FLD_START,ISEND_FLD_END

REAL(KIND=JPRB), TARGET :: PCOMBUFS_STACK(-1:YDCTX%ISENDCOUNT,MERGE (YDCTX%INSEND,0,NSTACK_MEMORY_TR/=0))
REAL(KIND=JPRB), TARGET :: PCOMBUFR_STACK(-1:YDCTX%IRECVCOUNT,MERGE (YDCTX%INRECV,0,NSTACK_MEMORY_TR/=0))

REAL(KIND=JPRB), ALLOCATABLE, TARGET, SAVE :: PCOMBUFS_HEAP(:,:)
REAL(KIND=JPRB), ALLOCATABLE, TARGET, SAVE :: PCOMBUFR_HEAP(:,:)

REAL(KIND=JPRB), POINTER, CONTIGUOUS :: PCOMBUFS(:,:)
REAL(KIND=JPRB), POINTER, CONTIGUOUS :: PCOMBUFR(:,:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE_BAR

!     ------------------------------------------------------------------

!*       0.    Some initializations
!              --------------------
ASSOCIATE (KSENDCOUNT=>YDCTX%ISENDCOUNT, KRECVCOUNT=>YDCTX%IRECVCOUNT, KNSEND=>YDCTX%INSEND, KNRECV=>YDCTX%INRECV, &
& KSENDTOT=>YDCTX%ISENDTOT, KRECVTOT=>YDCTX%IRECVTOT, KSEND=>YDCTX%ISEND, KRECV=>YDCTX%IRECV, KINDEX=>YDCTX%IINDEX, &
& KNDOFF=>YDCTX%INDOFF, KGPTRSEND =>YDCTX%IGPTRSEND, KSETAL=>YDCTX%ISETAL, KSETBL=>YDCTX%ISETBL, KSETWL=>YDCTX%ISETWL, &
& KSETVL=>YDCTX%ISETVL)

IF (NSTACK_MEMORY_TR == 0) THEN
  CALL ALLOCATE_HEAP_BUFFER(PCOMBUFS_HEAP, YDCTX%ISENDCOUNT, YDCTX%INSEND)
  CALL ALLOCATE_HEAP_BUFFER(PCOMBUFR_HEAP, YDCTX%IRECVCOUNT, YDCTX%INRECV)
    
! Now, force the OS to allocate this shared array right now, not when it starts to be used which is
! an OPEN-MP loop, that would cause a threads synchronization lock :
  IF (YDCTX%INSEND > 0 .AND. YDCTX%ISENDCOUNT >=-1) PCOMBUFS_HEAP(-1,1)=HUGE(1._JPRB)
  PCOMBUFS (-1:,1:) => PCOMBUFS_HEAP
  PCOMBUFR (-1:,1:) => PCOMBUFR_HEAP
ELSE
  PCOMBUFS (-1:,1:) => PCOMBUFS_STACK
  PCOMBUFR (-1:,1:) => PCOMBUFR_STACK
ENDIF

ITAG   = MTAGLG

IF (LHOOK) CALL DR_HOOK('TRLTOG_BAR',0,ZHOOK_HANDLE_BAR)
CALL GSTATS_BARRIER(762)
IF (LHOOK) CALL DR_HOOK('TRLTOG_BAR',1,ZHOOK_HANDLE_BAR)

CALL GSTATS(805,0)

IF (NTRANS_SYNC_LEVEL <= 0) THEN
  !...Receive loop.........................................................
  DO INR=1,KNRECV
    IRECV=KRECV(INR)
    CALL MPL_RECV(PCOMBUFR(-1:KRECVTOT(IRECV),INR), &
           & KSOURCE=NPRCIDS(IRECV), &
           & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQ_RECV(INR), &
           & KTAG=ITAG,CDSTRING='TRLTOG_COMM: NON-BLOCKING IRECV' )
  ENDDO
ENDIF

CALL GSTATS(805,1)

CALL GSTATS(1806,0)
LLINDER = .FALSE.
LLPGPONLY = .FALSE.
IF(PRESENT(KPTRGP))  LLINDER = .TRUE.
IF(PRESENT(PGP))     LLPGPONLY = .TRUE.
CALL ALLOCATE_TOG_VARS(YDTOGLVARS, KF_GP,NGPBLKS)
CALL INIT_TOG_VARS(YDTOGLVARS, KF_SCALARS_G, &
                   & PGP, PGPUV, PGP3A, PGP3B, PGP2)

CALL GSTATS(1806,1)


! Copy local contribution

IF(KSENDTOT(MYPROC) > 0 )THEN

  CALL INIT_TOG_OFF_VARS(YDCTX,YDTOGLVARS,KVSET,KPTRGP,KF_GP,NGPBLKS)

  CALL GSTATS(1604,0)

  CALL COPY_PGLAT(.TRUE.,YDCTX, YDTOGLVARS, NGPBLKS,IFLDS, LLPGPONLY, LLINDER,&                    
                    PGP=PGP, PGPUV=PGPUV, PGP3A=PGP3A, PGP3B=PGP3B, PGP2_RDLY=PGP2, &
                    & PGLAT_RDLY =PGLAT)

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
      PCOMBUFS((JFLD-ISEND_FLD_START)*ILEN+JL,INS) = PGLAT(JFLD,II)
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  PCOMBUFS(-1,INS) = 1
  PCOMBUFS(0,INS)  = KF_FS
  IF (NTRANS_SYNC_LEVEL <= 1) THEN
    CALL MPL_SEND(PCOMBUFS(-1:KSENDTOT(ISEND),INS),KDEST=NPRCIDS(ISEND),&
          & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQ_SEND(INS), &
          & KTAG=ITAG,CDSTRING='TRLTOG_COMM: NON-BLOCKING ISEND')
  ELSE
    CALL MPL_SEND(PCOMBUFS(-1:KSENDTOT(ISEND),INS),KDEST=NPRCIDS(ISEND),&
          & KMP_TYPE=JP_BLOCKING_BUFFERED, &
          & KTAG=ITAG,CDSTRING='TRLTOG_COMM: BLOCKING BUFFERED BSEND')
  ENDIF
ENDDO

!  Unpack loop.........................................................

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(INR,IRECV,ISETV,IFLD,JFLD,IPOS,JBLK,IFIRST,ILAST)
DO INR=1,KNRECV
  IRECV=KRECV(INR)
  
  ISETW(INR)=KSETWL(IRECV)
  ISETV=KSETVL(IRECV)
  IFLD = 0
  DO JFLD=1,KF_GP
    IF(KVSET(JFLD) == ISETV .OR. KVSET(JFLD) == -1 ) THEN
      IFLD = IFLD+1
      IFLDA(IFLD,INR)=JFLD
    ENDIF
  ENDDO
  IPOS = 0
  IPOSPLUS(INR)=0
  DO JBLK=1,NGPBLKS
    IFIRST = KGPTRSEND(1,JBLK,ISETW(INR))
    IF(IFIRST > 0) THEN
      ILAST = KGPTRSEND(2,JBLK,ISETW(INR))
      IJPOS(JBLK,INR)=IPOS
      IPOSPLUS(INR)=IPOSPLUS(INR)+(ILAST-IFIRST+1)
      IPOS=IPOS+(ILAST-IFIRST+1)
    ENDIF
  ENDDO
ENDDO
!$OMP END PARALLEL DO

DO JNR=1,KNRECV
  
  IF (NTRANS_SYNC_LEVEL <= 0) THEN
    CALL MPL_WAITANY(KREQUEST=IREQ_RECV(1:KNRECV),KINDEX=INR,&
          & CDSTRING='TRLTOG_COMM: WAIT FOR ANY RECEIVES')
  ELSE
    INR = JNR
    IRECV=KRECV(INR)
    CALL MPL_RECV(PCOMBUFR(-1:KRECVTOT(IRECV),INR), &
          & KSOURCE=NPRCIDS(IRECV), &
          & KMP_TYPE=JP_BLOCKING_STANDARD, &
          & KTAG=ITAG,CDSTRING='TRLTOG_COMM: BLOCKING RECV' )
  ENDIF

  CALL COPY_ZCOMBUF(.TRUE., YDCTX, YDTOGLVARS, INR, &
                      &  PCOMBUFR,NGPBLKS, LLPGPONLY, LLINDER,&
                      & IPOSPLUS, IFLDA,ISETW,IJPOS, &
                      & KPTRGP, &
                      PGP=PGP, PGPUV=PGPUV, PGP3A=PGP3A, PGP3B=PGP3B, PGP2=PGP2) 
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
