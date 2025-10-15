! (C) Copyright 1995- ECMWF.
! (C) Copyright 1995- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TRGTOL_MOD

PUBLIC TRGTOL
PRIVATE TRGTOL_PROLOG, TRGTOL_COMM

CONTAINS

SUBROUTINE TRGTOL(PGLAT,KF_FS,KF_GP,KF_SCALARS_G,KVSET,KPTRGP,&
 &PGP,PGPUV,PGP3A,PGP3B,PGP2)

!**** *TRGTOL * - head routine for transposition of grid point data from column
!                 structure to latitudinal. Reorganize data between
!                 grid point calculations and direct Fourier Transform

!**   Interface.
!     ----------
!        *call* *trgtol_prolog(...)

!        Explicit arguments :
!        --------------------

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
!        Original  : 18-Aug-2014 from trgtol
!        R. El Khatib 09-Sep-2020 NSTACK_MEMORY_TR
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DISTR       ,ONLY : D, NPRTRNS, NPROC
USE TPM_TRANS       ,ONLY : NGPBLKS
USE TRXTOY_CTX_MOD, ONLY : TRXTOY_CTX, TRXTOY_VARS

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KF_FS,KF_GP
REAL(KIND=JPRB),INTENT(OUT)   :: PGLAT(KF_FS,D%NLENGTF)
INTEGER(KIND=JPIM),INTENT(IN) :: KVSET(KF_GP)
INTEGER(KIND=JPIM),INTENT(IN) :: KF_SCALARS_G
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KPTRGP(:)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)     :: PGP(:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)     :: PGPUV(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)     :: PGP3A(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)     :: PGP3B(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)     :: PGP2(:,:,:)

TYPE (TRXTOY_CTX) :: YLCTX

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TRGTOL',0,ZHOOK_HANDLE)

CALL TRGTOL_PROLOG(KF_FS,KF_GP,KVSET,YLCTX)

CALL TRGTOL_COMM(PGLAT,KF_FS,KF_GP,KF_SCALARS_G,KVSET, &
                 & KPTRGP,PGP,PGPUV,PGP3A,PGP3B,PGP2,YLCTX)

IF (LHOOK) CALL DR_HOOK('TRGTOL',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE TRGTOL

SUBROUTINE TRGTOL_PROLOG(KF_FS,KF_GP,KVSET,YDCTX)

!**** *TRGTOL_PROLOG * - prolog for transposition of grid point data from column
!                 structure to latitudinal. Reorganize data between
!                 grid point calculations and direct Fourier Transform
!                 the purpose is essentially 
!                 to compute the size of communication buffers in order to enable
!                 the use of automatic arrays later.


!     Purpose.
!     --------


!**   Interface.
!     ----------
!        *call* *trgtol_prolog(...)

!        Explicit arguments :
!        --------------------

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
!        Original  : 18-Aug-2014 from trgtol
!     ------------------------------------------------------------------



USE PARKIND1  ,ONLY : JPIM 

USE TPM_DISTR       ,ONLY : D, NPRTRNS, MYSETW, MYPROC, NPROC
USE TPM_TRANS       ,ONLY : NGPBLKS

USE INIGPTR_MOD     ,ONLY : INIGPTR
USE PE2SET_MOD      ,ONLY : PE2SET
USE TRXTOY_CTX_MOD  ,ONLY : TRXTOY_CTX,ALLOCATE_CTX

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KF_FS,KF_GP
TYPE (TRXTOY_CTX), INTENT(INOUT) :: YDCTX
INTEGER(KIND=JPIM),INTENT(IN) :: KVSET(KF_GP)

INTEGER(KIND=JPIM) :: IGPTRRECV(NPRTRNS)
INTEGER(KIND=JPIM) :: IFIRSTLAT, IGL, IGLL, ILASTLAT, JROC, IPOS, ISETB, ISETA
INTEGER(KIND=JPIM) :: ISETV, J, JFLD, JGL, JL, ISETW,  INDOFFX

!     ------------------------------------------------------------------
!*       0.    Some initializations
!              --------------------

CALL ALLOCATE_CTX(YDCTX, NPROC,D%NLENGTF,NGPBLKS,NPRTRNS)

CALL GSTATS(1805,0)

CALL INIGPTR(YDCTX%IGPTRSEND,IGPTRRECV)

INDOFFX  = 0
YDCTX%INRECV = 0
YDCTX%INSEND = 0

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
  YDCTX%ISENDTOT(JROC) = IGPTRRECV(ISETW)*IPOS

  IF( JROC /= MYPROC) THEN
    IF(YDCTX%ISENDTOT(JROC) > 0) THEN
      YDCTX%INSEND = YDCTX%INSEND+1
      YDCTX%ISEND(YDCTX%INSEND)=JROC
    ENDIF
  ENDIF

  IFIRSTLAT = MAX(D%NPTRLS(MYSETW),D%NFRSTLAT(ISETA))
  ILASTLAT  = MIN(D%NPTRLS(MYSETW)+D%NULTPP(MYSETW)-1,D%NLSTLAT(ISETA))

  IPOS = 0
  DO JGL=IFIRSTLAT,ILASTLAT
    IGL  = D%NPTRFRSTLAT(ISETA)+JGL-D%NFRSTLAT(ISETA)
    IPOS = IPOS+D%NONL(IGL,ISETB)
  ENDDO

  YDCTX%IRECVTOT(JROC) = IPOS*KF_FS
  IF(YDCTX%IRECVTOT(JROC) > 0 .AND. MYPROC /= JROC) THEN
    YDCTX%INRECV = YDCTX%INRECV + 1
    YDCTX%IRECV(YDCTX%INRECV)=JROC
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

CALL GSTATS(1805,1)

END SUBROUTINE TRGTOL_PROLOG

SUBROUTINE TRGTOL_COMM(PGLAT,KF_FS,KF_GP,KF_SCALARS_G,KVSET,&
 & KPTRGP,PGP,PGPUV,PGP3A,PGP3B,PGP2,YDCTX)

!**** *TRGTOL_COMM * - transposition of grid point data from column
!                 structure to latitudinal. Reorganize data between
!                 grid point calculations and direct Fourier Transform

!     Purpose.
!     --------


!**   Interface.
!     ----------
!        *call* *trgtol(...)

!        Explicit arguments :
!        --------------------
!           PGLAT    -  Latitudinal data ready for direct FFT (output)
!           PGP    -  Blocked grid point data    (input)

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
!        Original: 95-10-01
!        D.Dent  : 97-08-04   Reorganisation to allow
!                             NPRTRV to differ from NPRGPEW
!                : 98-06-17   add mailbox control logic (from TRLTOM)
!        =99-03-29= Mats Hamrud and Deborah Salmond
!                   JUMP in FFT's changed to 1
!                   KINDEX introduced and PCOMBUF not used for same PE
!         01-11-23  Deborah Salmond and John Hague
!                    LIMP_NOOLAP Option for non-overlapping message passing
!                    and buffer packing
!         01-12-18  Peter Towers
!                   Improved vector performance of GTOL_PACK,GTOL_UNPACK
!         03-04-02  G. Radnoti: call barrier always when nproc>1
!         08-01-01  G.Mozdzynski: cleanup
!         09-01-02  G.Mozdzynski: use non-blocking recv and send
!        R. El Khatib 09-Sep-2020 64 bits addressing for PGLAT
!     ------------------------------------------------------------------



USE PARKIND1  ,ONLY : JPIM     ,JPRB    ,JPIB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE MPL_MODULE  ,ONLY : MPL_RECV, MPL_SEND, MPL_WAIT, JP_NON_BLOCKING_STANDARD, MPL_WAITANY, &
     & JP_BLOCKING_STANDARD, MPL_BARRIER, JP_BLOCKING_BUFFERED

USE TPM_GEN         ,ONLY : NOUT, NTRANS_SYNC_LEVEL, NSTACK_MEMORY_TR
USE TPM_DISTR       ,ONLY : D, NPRCIDS, NPRTRNS, MTAGGL,  &
     &                      MYSETV, MYSETW, MYPROC, NPROC
USE TPM_TRANS       ,ONLY : LDIVGP, LGPNORM, LSCDERS, LUVDER, LVORGP, NGPBLKS

USE PE2SET_MOD      ,ONLY : PE2SET
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
USE TRXTOY_CTX_MOD  ,ONLY : TRXTOY_CTX, TRXTOY_VARS, &
                             &ALLOCATE_TOG_VARS, ALLOCATE_HEAP_BUFFER,&
                             &INIT_TOG_VARS, INIT_TOG_OFF_VARS, &
                             &COPY_ZCOMBUF, COPY_PGLAT


IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KF_FS,KF_GP
REAL(KIND=JPRB),INTENT(OUT)   :: PGLAT(KF_FS,D%NLENGTF)
INTEGER(KIND=JPIM),INTENT(IN) :: KVSET(KF_GP)
INTEGER(KIND=JPIM),INTENT(IN) :: KF_SCALARS_G
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KPTRGP(:)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)     :: PGP(:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)     :: PGPUV(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)     :: PGP3A(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)     :: PGP3B(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)     :: PGP2(:,:,:)

TYPE(TRXTOY_CTX), INTENT(INOUT), TARGET :: YDCTX
! LOCAL VARIABLES
TYPE(TRXTOY_VARS) :: YDTOGLVARS
INTEGER(KIND=JPIM) :: IPOSPLUS(YDCTX%INSEND)
INTEGER(KIND=JPIM) :: ISETW(YDCTX%INSEND)
INTEGER(KIND=JPIM) :: IJPOS(NGPBLKS,YDCTX%INSEND)
INTEGER(KIND=JPIM) :: IFLDA(KF_GP,YDCTX%INSEND)
INTEGER(KIND=JPIM) :: IREQ_SEND(NPROC)
INTEGER(KIND=JPIM) :: IREQ_RECV(NPROC)

!     LOCAL LOGICAL SCALARS
LOGICAL   :: LLPGPONLY, LLINDER
LOGICAL   :: LLUV(KF_GP),LLGP2(KF_GP),LLGP3A(KF_GP),LLGP3B(KF_GP)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: IFIRST, ILAST, ILEN, IPOS, ISETA, ISETB, IRECV, ISETV
INTEGER(KIND=JPIM) :: ISEND, ITAG, JBLK, JFLD, JK, JL, IFLD, II, IFLDS, INS, INR
INTEGER(KIND=JPIM) :: JJ,JI,IFLDT, J

INTEGER(KIND=JPIB) :: JFLD64

INTEGER(KIND=JPIM) :: IUVLEVS(KF_GP),IUVPARS(KF_GP),IGP2PARS(KF_GP)
INTEGER(KIND=JPIM) :: IGP3APARS(KF_GP),IGP3ALEVS(KF_GP),IGP3BPARS(KF_GP),IGP3BLEVS(KF_GP)
INTEGER(KIND=JPIM) :: IUVPAR,IUVLEV,IGP2PAR,IGP3ALEV,IGP3APAR,IGP3BLEV,IGP3BPAR,IPAROFF
INTEGER(KIND=JPIM) :: IOFF,IOFF1,IOFFNS,IOFFEW,J1,J2, JNR
INTEGER(KIND=JPIM) :: IFLDOFF(KF_FS)
INTEGER(KIND=JPIM) :: ISEND_FLD_START,ISEND_FLD_END
INTEGER(KIND=JPIM) :: IRECV_FLD_START,IRECV_FLD_END
INTEGER(KIND=JPIM) :: IGPTROFF(NGPBLKS)
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
& KNDOFF=>YDCTX%INDOFF, KGPTRSEND =>YDCTX%IGPTRSEND, KSETAL=>YDCTX%ISETAL, KSETBL=>YDCTX%ISETBL, KSETWL=>YDCTX%ISETWL, &
& KSETVL=>YDCTX%ISETVL)

IF (NSTACK_MEMORY_TR == 0) THEN
  CALL ALLOCATE_HEAP_BUFFER(ZCOMBUFS_HEAP, YDCTX%ISENDCOUNT, YDCTX%INSEND)
  CALL ALLOCATE_HEAP_BUFFER(ZCOMBUFR_HEAP, YDCTX%IRECVCOUNT, YDCTX%INRECV)

! Now, force the OS to allocate this shared array right now, not when it starts to be used which is
! an OPEN-MP loop, that would cause a threads synchronization lock :
  IF (YDCTX%INSEND > 0 .AND. YDCTX%ISENDCOUNT >=-1) ZCOMBUFS_HEAP(-1,1)=HUGE(1._JPRB)
  ZCOMBUFS (-1:,1:) => ZCOMBUFS_HEAP
  ZCOMBUFR (-1:,1:) => ZCOMBUFR_HEAP
ELSE
  ZCOMBUFS (-1:,1:) => ZCOMBUFS_STACK
  ZCOMBUFR (-1:,1:) => ZCOMBUFR_STACK
ENDIF

ITAG = MTAGGL

IF (LHOOK) CALL DR_HOOK('TRGTOL_BAR',0,ZHOOK_HANDLE_BAR)
CALL GSTATS_BARRIER(761)
IF (LHOOK) CALL DR_HOOK('TRGTOL_BAR',1,ZHOOK_HANDLE_BAR)

IF(.NOT.LGPNORM)THEN
  CALL GSTATS(803,0)
ELSE
  CALL GSTATS(804,0)
ENDIF

IF (NTRANS_SYNC_LEVEL <= 0) THEN
  !  Receive loop.........................................................
  DO INR=1,KNRECV
    IRECV=KRECV(INR)
    CALL MPL_RECV(ZCOMBUFR(-1:KRECVTOT(IRECV),INR), &
        & KSOURCE=NPRCIDS(IRECV), &
        & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQ_RECV(INR), &
        & KTAG=ITAG,CDSTRING='TRGTOL_COMM: NON-BLOCKING IRECV' )
  ENDDO
ENDIF

IF(.NOT.LGPNORM)THEN
  CALL GSTATS(803,1)
ELSE
  CALL GSTATS(804,1)
ENDIF

CALL GSTATS(1805,0)
LLINDER = .FALSE.
LLPGPONLY = .FALSE.
IF(PRESENT(KPTRGP))  LLINDER = .TRUE.
IF(PRESENT(PGP))     LLPGPONLY = .TRUE.
CALL ALLOCATE_TOG_VARS(YDTOGLVARS, KF_GP,NGPBLKS)
CALL INIT_TOG_VARS(YDTOGLVARS, KF_SCALARS_G, &
                   & PGP, PGPUV, PGP3A, PGP3B, PGP2)

CALL GSTATS(1805,1)


! Copy local contribution

IF(KSENDTOT(MYPROC) > 0 )THEN

  CALL INIT_TOG_OFF_VARS(YDCTX,YDTOGLVARS,KVSET,KPTRGP,KF_GP,NGPBLKS)
  
  CALL GSTATS(1601,0)
  
  CALL COPY_PGLAT(.FALSE.,YDCTX, YDTOGLVARS, NGPBLKS,IFLDS, LLPGPONLY, LLINDER,&                    
                    & PGP_RDLY=PGP, PGPUV_RDLY=PGPUV, PGP3A_RDLY=PGP3A, PGP3B_RDLY=PGP3B, PGP2_RDLY=PGP2, &
                    & PGLAT =PGLAT)

  CALL GSTATS(1601,1)

ENDIF


! Now overlapping buffer packing/unpacking with sends/waits
! Time as if all communications to avoid double accounting

IF(.NOT.LGPNORM)THEN
  CALL GSTATS(803,0)
ELSE
  CALL GSTATS(804,0)
ENDIF

!....Pack+send loop.........................................................

!$OMP PARALLEL PRIVATE(JBLK,IFIRST,ILAST,ISEND_FLD_START,ISEND_FLD_END,INS,ISEND,ISETA,ISETB,&
!$OMP&                 ISETV,IFLD,IFLDT,IPOS,JFLD,JK,JJ,JI)
!$OMP DO SCHEDULE(STATIC)
DO INS=1,KNSEND
  ISEND=KSEND(INS)
  CALL PE2SET(ISEND,ISETA,ISETB,ISETW(INS),ISETV)
  IFLD = 0
  IPOS = 0
  IPOSPLUS(INS)=0
  DO JFLD=1,KF_GP
    IF(KVSET(JFLD) == ISETV .OR. KVSET(JFLD) == -1 ) THEN
      IFLD = IFLD+1
      IFLDA(IFLD,INS)=JFLD
    ENDIF
  ENDDO

  DO JBLK=1,NGPBLKS
    IFIRST = KGPTRSEND(1,JBLK,ISETW(INS))
    IF(IFIRST > 0) THEN
      ILAST = KGPTRSEND(2,JBLK,ISETW(INS))
      IJPOS(JBLK,INS)=IPOS
      IPOSPLUS(INS)=IPOSPLUS(INS)+(ILAST-IFIRST+1)
      IPOS=IPOS+(ILAST-IFIRST+1)
    ENDIF
  ENDDO

  ZCOMBUFS(-1,INS) = 1
  ZCOMBUFS(0,INS) = IFLD

ENDDO
!$OMP END DO
DO INS=1,KNSEND
  
  CALL COPY_ZCOMBUF(.FALSE., YDCTX, YDTOGLVARS, INS, &
                      &  ZCOMBUFS,NGPBLKS, LLPGPONLY, LLINDER,&
                      & IPOSPLUS, IFLDA,ISETW,IJPOS, &
                      & KPTRGP, &
                      PGP_RDLY=PGP, PGPUV_RDLY=PGPUV, PGP3A_RDLY=PGP3A, PGP3B_RDLY=PGP3B, PGP2_RDLY=PGP2)
  
ENDDO
!$OMP END PARALLEL

DO INS=1,KNSEND
  ISEND=KSEND(INS)
  IF (NTRANS_SYNC_LEVEL <= 1) THEN
    CALL MPL_SEND(ZCOMBUFS(-1:YDCTX%ISENDTOT(ISEND),INS),KDEST=NPRCIDS(ISEND), &
          & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQ_SEND(INS), &
          & KTAG=ITAG,CDSTRING='TRGTOL_COMM: NON-BLOCKING ISEND' )
  ELSE
    CALL MPL_SEND(ZCOMBUFS(-1:YDCTX%ISENDTOT(ISEND),INS),KDEST=NPRCIDS(ISEND), &
          & KMP_TYPE=JP_BLOCKING_BUFFERED, &
          & KTAG=ITAG,CDSTRING='TRGTOL_COMM: BLOCKING BUFFERED BSEND' )
  ENDIF

ENDDO

!  Unpack loop.........................................................

DO JNR=1,KNRECV

  IF (NTRANS_SYNC_LEVEL <= 0) THEN
    CALL MPL_WAITANY(KREQUEST=IREQ_RECV(1:KNRECV),KINDEX=INR,&
          & CDSTRING='TRGTOL_COMM: WAIT FOR ANY RECEIVES')
  ELSE
    INR = JNR
    IRECV=KRECV(INR)
    CALL MPL_RECV(ZCOMBUFR(-1:YDCTX%IRECVTOT(IRECV),INR), &
          & KSOURCE=NPRCIDS(IRECV), &
          & KMP_TYPE=JP_BLOCKING_STANDARD, &
          & KTAG=ITAG,CDSTRING='TRGTOL_COMM: BLOCKING RECV' )
  ENDIF

  IRECV=KRECV(INR)
  ILEN = KRECVTOT(IRECV)/KF_FS
  IRECV_FLD_START = ZCOMBUFR(-1,INR)
  IRECV_FLD_END   = ZCOMBUFR(0,INR)
  !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JL,II,JFLD)
  DO JL=1,ILEN
    II = KINDEX(KNDOFF(IRECV)+JL)
    DO JFLD=IRECV_FLD_START,IRECV_FLD_END
      PGLAT(JFLD,II) = ZCOMBUFR(JL+(JFLD-IRECV_FLD_START)*ILEN,INR)
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  IPOS = ILEN*(IRECV_FLD_END-IRECV_FLD_START+1)
ENDDO

IF (NTRANS_SYNC_LEVEL <= 1) THEN
  IF(KNSEND > 0) THEN
    CALL MPL_WAIT(KREQUEST=IREQ_SEND(1:KNSEND),CDSTRING='TRGTOL_COMM: WAIT FOR ISENDS')
  ENDIF
ENDIF

IF (NTRANS_SYNC_LEVEL >= 1) THEN
  CALL MPL_BARRIER(CDSTRING='TRGTOL_COMM: BARRIER AT END')
ENDIF

IF(.NOT.LGPNORM)THEN
  CALL GSTATS(803,1)
ELSE
  CALL GSTATS(804,1)
ENDIF

CALL GSTATS_BARRIER2(761)

END ASSOCIATE

END SUBROUTINE TRGTOL_COMM

END MODULE TRGTOL_MOD
