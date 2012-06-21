MODULE TRLTOM_MOD
CONTAINS
SUBROUTINE TRLTOM(PFBUF_IN,PFBUF,KFIELD)

!**** *TRLTOM * - transposition in Fourierspace

!     Purpose.
!     --------
!              Transpose Fourier coefficients from partitioning
!              over latitudes to partitioning over wave numbers
!              This is done between inverse Legendre Transform
!              and inverse FFT.
!              This is the inverse routine of TRMTOL.

!**   Interface.
!     ----------
!        *CALL* *TRLTOM(...)*

!        Explicit arguments : PFBUF  - Fourier coefficient buffer. It is
!        --------------------          used for both input and output.

!                             KFIELD - Number of fields communicated

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
!        Original : 95-10-01
!        Modified : 97-06-18 G. Mozdzynski - control MPI mailbox use
!                                            (NCOMBFLEN) for nphase.eq.1
!        Modified : 99-05-28  D.Salmond - Optimise copies. 
!        Modified : 00-02-02  M.Hamrud  - Remove NPHASE
!        D.Salmond : 01-11-23 LIMP_NOOLAP Option for non-overlapping message
!                             passing and buffer packing
!        Y.Seity   : 07-08-30 Add barrier synchonisation under LSYNC_TRANS
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE MPL_MODULE

USE TPM_DISTR
USE TPM_GEN

USE SET2PE_MOD
USE MYSENDSET_MOD
USE MYRECVSET_MOD
USE ABORT_TRANS_MOD

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELD
REAL(KIND=JPRB)   ,INTENT(INOUT)  :: PFBUF(:)
REAL(KIND=JPRB)   ,INTENT(INOUT)  :: PFBUF_IN(:)

INTEGER(KIND=JPIM) :: ISENT    (NPRTRW-1)
INTEGER(KIND=JPIM) :: IRCVD    (NPRTRW-1)
INTEGER(KIND=JPIM) :: ISENDTOT (NPRTRW-1)
INTEGER(KIND=JPIM) :: IRECVTOT (NPRTRW-1)
INTEGER(KIND=JPIM) :: ISENDREQ (NPRTRW-1)
INTEGER(KIND=JPIM) :: IRECVREQ (NPRTRW-1)
INTEGER(KIND=JPIM) :: ILENS(NPRTRW),IOFFS(NPRTRW),ILENR(NPRTRW),IOFFR(NPRTRW)

INTEGER(KIND=JPIM) :: ICOMBFLENP, ILEN, ILREC, &
             &IRECV, IRECVSET, ISEND,ISTP, &
             &ISENDSET, ISTA, ITAG, J, &
             &INUMSENT

LOGICAL :: LLDONE
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------


!     TRANSPOSITION
!     -------------

! Set maximum transfer length

IF (LHOOK) CALL DR_HOOK('TRLTOM',0,ZHOOK_HANDLE)

  CALL GSTATS_BARRIER(763)

IF( NPRTRW > 1 )THEN
  ICOMBFLENP = NCOMBFLEN/(3*NPRTRW-4)
ENDIF

! Mark all processors as having no outstanding receives,
! nothing sent or received.

DO J=1,NPRTRW-1
  ISENDSET = MYSENDSET(NPRTRW,MYSETW,J)
  IRECVSET = MYRECVSET(NPRTRW,MYSETW,J)
  ISENT(J) = 0
  IRCVD(J) = 0
  ISENDTOT(J) = D%NLTSGTB(ISENDSET)*KFIELD
  IRECVTOT(J) = D%NLTSFTB(IRECVSET)*KFIELD
ENDDO
LLDONE = .FALSE.
ITAG = MTAGLM


! Loop over the number of processors we need to communicate with

IF(LIMP_NOOLAP)THEN

! Send/Recv loop.............................................................

  DO J=1,NPRTRW
    ILENS(J) = D%NLTSGTB(J)*KFIELD
    IOFFS(J) = D%NSTAGT1B(D%MSTABF(J))*KFIELD
    ILENR(J) = D%NLTSFTB(J)*KFIELD
    IOFFR(J) = D%NSTAGT1B(J)*KFIELD
  ENDDO
  CALL GSTATS(806,0)

!ajout mpl_barrier
  IF (LSYNC_TRANS) THEN
    CALL MPL_BARRIER(CDSTRING='TRLTOM:')
  ENDIF

  CALL MPL_ALLTOALLV(PSENDBUF=PFBUF_IN,KSENDCOUNTS=ILENS,&
   & PRECVBUF=PFBUF,KRECVCOUNTS=ILENR,KSENDDISPL=IOFFS,KRECVDISPL=IOFFR,&
   & KCOMM=MPL_ALL_MS_COMM,CDSTRING='TRLTOM:')
  
  CALL GSTATS(806,1)

ELSE
! Copy local part of PFBUF_IN to PFBUF

  ILEN = D%NLTSGTB(MYSETW)*KFIELD
  ISTA = D%NSTAGT1B(MYSETW)*KFIELD+1
  CALL GSTATS(1607,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(J)
  DO J=ISTA,ISTA+ILEN-1
    PFBUF(J) = PFBUF_IN(J)
  ENDDO
!$OMP END PARALLEL DO
  CALL GSTATS(1607,1)

  DO WHILE( .NOT. LLDONE )
  
    LLDONE = .TRUE.
  
  
! Send loop.............................................................
  
!  For immediate send/recv, post receives first
  
    IF(LIMP) THEN
      DO J=1,NPRTRW-1
        IF( IRECVTOT(J)-IRCVD(J) > 0 )THEN
          IRECVSET = MYRECVSET(NPRTRW,MYSETW,J)
          CALL SET2PE(IRECV,0,0,IRECVSET,MYSETV)
          ISTA = D%NSTAGT1B(IRECVSET)*KFIELD+IRCVD(J)+1
          ILEN = MIN(ICOMBFLENP,IRECVTOT(J)-IRCVD(J))
          ISTP = ISTA+ILEN-1
          CALL MPL_RECV(PFBUF(ISTA:ISTP),KSOURCE=NPRCIDS(IRECV),KTAG=ITAG, &
          &    KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IRECVREQ(J), &
          &    CDSTRING='TRLTOM:' )
        ENDIF
      ENDDO
    ENDIF
  
    INUMSENT = 0
    DO J=1,NPRTRW-1
  
! Check if there is no outstanding receive, and we have data to send
  
      IF( ISENDTOT(J)-ISENT(J) > 0 )THEN
  
        LLDONE = .FALSE.
        ISENDSET = MYSENDSET(NPRTRW,MYSETW,J)
        CALL SET2PE(ISEND,0,0,ISENDSET,MYSETV)
        ILEN = MIN(ICOMBFLENP,ISENDTOT(J)-ISENT(J))
        ISTA = D%NSTAGT1B(D%MSTABF(ISENDSET))*KFIELD+ISENT(J)+1
        ISTP = ISTA+ILEN-1
  
        IF(LIMP) THEN
          INUMSENT = INUMSENT+1
          CALL MPL_SEND(PFBUF_IN(ISTA:ISTP),KDEST=NPRCIDS(ISEND),&
      &    KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(INUMSENT), &
           &KTAG=ITAG,CDSTRING='TRLTOM:')
        ELSE
          CALL MPL_SEND(PFBUF_IN(ISTA:ISTP),KDEST=NPRCIDS(ISEND),&
           &KTAG=ITAG,CDSTRING='TRLTOM:')
        ENDIF
  
        ISENT(J) = ISENT(J)+ILEN
        
      ENDIF
  
    ENDDO
  
  
!  Receive loop.........................................................
  
  
    DO J=1,NPRTRW-1
      IF( IRECVTOT(J)-IRCVD(J) > 0 )THEN
        LLDONE = .FALSE.
        ILEN = MIN(ICOMBFLENP,IRECVTOT(J)-IRCVD(J))
        IF( .NOT. LIMP ) THEN
          IRECVSET = MYRECVSET(NPRTRW,MYSETW,J)
          CALL SET2PE(IRECV,0,0,IRECVSET,MYSETV)
          ISTA = D%NSTAGT1B(IRECVSET)*KFIELD+IRCVD(J)+1
          ISTP = ISTA+ILEN-1
          CALL MPL_RECV(PFBUF(ISTA:ISTP),KSOURCE=NPRCIDS(IRECV),&
           &KTAG=ITAG,KOUNT=ILREC,CDSTRING='TRLTOM:')
        ELSE
!       For LIMP=true, find message length
          CALL MPL_WAIT(PFBUF,KREQUEST=IRECVREQ(J), &
           & KOUNT=ILREC,CDSTRING='TRLTOM: LIMP WAITS ' )
        ENDIF
        IF( ILREC /= ILEN )THEN
          CALL ABORT_TRANS('TRLTOM:RECEIVED MESSAGE LENGTH ERROR')
        ENDIF
        IRCVD(J) = IRCVD(J)+ILEN
      ENDIF
    ENDDO
  
!       For LIMP=true, wait for sends to complete
    IF(LIMP) THEN
      IF( INUMSENT > 0 ) THEN
        CALL MPL_WAIT(PFBUF,KREQUEST=ISENDREQ(1:INUMSENT), &
         & CDSTRING='TRLTOM: ERROR IN MPL_WAIT FOR SENDS')
      ENDIF
    ENDIF
  
  ENDDO
! Perform barrier synchronisation to guarantee all processors have
! completed communication

  IF( NPRTRW > 1 )THEN
    CALL MPL_BARRIER(CDSTRING='TRLTOM:')
  ENDIF
ENDIF

IF (LHOOK) CALL DR_HOOK('TRLTOM',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE TRLTOM
END MODULE TRLTOM_MOD
