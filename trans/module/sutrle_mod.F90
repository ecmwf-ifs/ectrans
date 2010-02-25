MODULE SUTRLE_MOD
CONTAINS
SUBROUTINE SUTRLE(PNM,KGL,KLOOP,KEXPECT)

!**** *sutrle * - transposition of Legendre polynomials during set-up

!     Purpose.
!     --------
!           transposition of Legendre polynomials during set-up

!**   Interface.
!     ----------
!        *call* *sutrle(pnm)

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
!        MPP Group *ECMWF*

!     Modifications.
!     --------------
!        Original : 95-10-01
!        P.Towers : 10-01-12  Corrected over allocation of ZSNDBUF (XT4 fix)
!     ------------------------------------------------------------------


USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE MPL_MODULE

USE TPM_GEN
USE TPM_DIM
USE TPM_DISTR
USE TPM_FIELDS
USE SET2PE_MOD
USE ABORT_TRANS_MOD

IMPLICIT NONE

REAL(KIND=JPRB),INTENT(IN) :: PNM(:)
INTEGER(KIND=JPIM),INTENT(IN) :: KGL
INTEGER(KIND=JPIM),INTENT(IN) :: KLOOP
INTEGER(KIND=JPIM),INTENT(INOUT) :: KEXPECT(:)
!     LOCAL 

REAL(KIND=JPRB), ALLOCATABLE :: ZSNDBUF(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZRCVBUF(:)
INTEGER(KIND=JPIM) :: ILREC, IM, IPOS, &
             & IRECVSET, IRECV, ISEND, ISENDSET, ITAG,ISENDSIZE, &
             & JM, JMLOC, JN, JROC ,IOFFT, IOFFG,IGL,ISREQ
INTEGER(KIND=JPIM) :: ISENDREQ(NPRTRW)

!     ------------------------------------------------------------------

!*       0.    Some initializations.
!              ---------------------

! Perform barrier synchronisation to guarantee all processors have
! completed all previous communication

IF( NPROC > 1 .AND. KLOOP ==1)THEN
  CALL GSTATS(783,0)
  CALL MPL_BARRIER(CDSTRING='SUTRLE:')
  CALL GSTATS(783,1)
ENDIF

!*     Calculate send buffer size

ISENDSIZE=0
DO JROC=1,NPRTRW-1

  ISEND = MYSETW-JROC
  IF (ISEND <= 0)     ISEND = ISEND+NPRTRW
  ISENDSET = ISEND
  CALL SET2PE(ISEND,0,0,ISENDSET,MYSETV)

  IF(KGL > 0) THEN
    IPOS = 1
    DO JM=0,R%NSMAX
      IF (ISENDSET == D%NPROCM(JM) ) IPOS = IPOS + R%NTMAX-JM+2 
    ENDDO
    ISENDSIZE = MAX(IPOS,ISENDSIZE)
  ENDIF
ENDDO

ALLOCATE (ZSNDBUF(ISENDSIZE,NPRTRW))
ALLOCATE (ZRCVBUF(R%NSPOLEG+1))
!*   copy data to be sent into zsndbuf

ISREQ = 0
DO JROC=1,NPRTRW-1

  ITAG = MTAGLETR+KLOOP

!*     Define PE to which data have to be sent and PE from which
!*     data have to be received

  CALL GSTATS(801,0)
  ISEND = MYSETW-JROC
  IRECV = MYSETW+JROC
  IF (ISEND <= 0)     ISEND = ISEND+NPRTRW
  IF (IRECV > NPRTRW) IRECV = IRECV-NPRTRW
  IRECVSET = IRECV
  ISENDSET = ISEND
  CALL SET2PE(ISEND,0,0,ISENDSET,MYSETV)
  CALL SET2PE(IRECV,0,0,IRECVSET,MYSETV)


!*   send message 
  
  
  IF(KGL > 0) THEN
    ZSNDBUF(1,ISENDSET) = kgl
    IPOS = 1
    DO JM=0,R%NSMAX
      IF (ISENDSET == D%NPROCM(JM) ) THEN
        DO JN=1,R%NTMAX-JM+2
          IPOS = IPOS + 1
          ZSNDBUF(IPOS,ISENDSET) = PNM(D%NPMG(JM)+JN)
        ENDDO
      ENDIF
    ENDDO
    ISENDSIZE = IPOS
    ISREQ = ISREQ+1
    CALL MPL_SEND(ZSNDBUF(1:ISENDSIZE,ISENDSET),KDEST=NPRCIDS(ISEND), &
     &KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(ISREQ),&
     & KTAG=ITAG,CDSTRING='SUTRLE:')
  ENDIF

  CALL GSTATS(801,1)

  ILREC = 0
  IF (D%NUMP > 0.AND. KEXPECT(IRECVSET) > 0) THEN

!*   receive message (if not empty)

    CALL GSTATS(801,0)
    CALL MPL_RECV(ZRCVBUF,KSOURCE=NPRCIDS(IRECV), &
     & KTAG=ITAG,KOUNT=ILREC,CDSTRING='SUTRLE:')

!*   copy data from buffer to f%rpnm
    
    IGL = ZRCVBUF(1)
    IPOS = 1
    IOFFT = D%NPMT(IGL) 
    DO JMLOC=1,D%NUMP
      IM = D%MYMS(JMLOC)
      IOFFT = D%NPMT(IM)
      DO JN=1,R%NTMAX-IM+2
        IPOS = IPOS + 1
        F%RPNM(IGL,IOFFT+JN) = ZRCVBUF(IPOS)
      ENDDO
    ENDDO
!*    check received message length

    IF (ILREC /= IPOS) THEN
      WRITE(NOUT,*)' SUTRLE: ILREC,IPOS,NCOMBLEN ',ILREC,IPOS,NCOMBFLEN
      CALL ABORT_TRANS(' SUTRLE:RECEIVED MESSAGE LENGTH DOES NOT MATCH')
    ENDIF
    KEXPECT(IRECVSET) = KEXPECT(IRECVSET)-1
    CALL GSTATS(801,1)

  ENDIF
ENDDO
IF(ISREQ > 0) THEN
  CALL MPL_WAIT(ZSNDBUF(:,1),KREQUEST=ISENDREQ(1:ISREQ), &
   & CDSTRING='SUTRLE: WAIT')
ENDIF
!*    copy data from pnm to rpnm

IF(KEXPECT(MYSETW) > 0) THEN
  CALL GSTATS(1803,0)
  DO JMLOC=1,D%NUMP
    IM = D%MYMS(JMLOC)
    IOFFT = D%NPMT(IM)
    IOFFG = D%NPMG(IM)
    DO JN=1,R%NTMAX-IM+2
      F%RPNM(KGL,IOFFT+JN) = PNM(IOFFG+JN)
    ENDDO
  ENDDO
  KEXPECT(MYSETW) = KEXPECT(MYSETW)-1
  CALL GSTATS(1803,1)
ENDIF
DEALLOCATE (ZSNDBUF)
DEALLOCATE (ZRCVBUF)

END SUBROUTINE SUTRLE
END MODULE SUTRLE_MOD
