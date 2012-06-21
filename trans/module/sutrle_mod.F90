MODULE SUTRLE_MOD
CONTAINS
SUBROUTINE SUTRLE(PNM)

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

REAL(KIND=JPRB),INTENT(IN) :: PNM(R%NSPOLEG,D%NLEI3D)

!     LOCAL 

REAL(KIND=JPRB), ALLOCATABLE :: ZCOMBUF(:)
REAL(KIND=JPRB), POINTER     :: ZPNM(:,:)
INTEGER(KIND=JPIM) :: IGLLOC, ILREC, IM, INENTR, IPOS, &
             &IRECSET, IRECV, ISEND, ISENDSET, ITAG, &
             &JGL, JGLLOC, JM, JMLOC, JN, JROC ,IOFFT, IOFFG

LOGICAL :: LLADMSG, LLEXACT

!     ------------------------------------------------------------------

!*       0.    Some initializations.
!              ---------------------
!! Workaround for obscure unwillingness to vectorize on VPP
ZPNM => F%RPNM

! Perform barrier synchronisation to guarantee all processors have
! completed all previous communication

IF( NPROC > 1 )THEN
  CALL GSTATS(783,0)
  CALL MPL_BARRIER(CDSTRING='SUTRLE:')
  CALL GSTATS(783,1)
ENDIF

ALLOCATE (ZCOMBUF(NCOMBFLEN))

DO JROC=1,NPRTRW-1

  LLADMSG = .FALSE.
  ITAG = MTAGLETR

!*     Define PE to which data have to be sent and PE from which
!*     data have to be received

  CALL GSTATS(801,0)
  ISEND = MYSETW-JROC
  IRECV = MYSETW+JROC
  IF (ISEND <= 0)     ISEND = ISEND+NPRTRW
  IF (IRECV > NPRTRW) IRECV = IRECV-NPRTRW
  IRECSET = IRECV
  ISENDSET = ISEND
  CALL SET2PE(ISEND,0,0,ISEND,MYSETV)
  CALL SET2PE(IRECV,0,0,IRECV,MYSETV)

!*   copy data to be sent into zcombuf

  IPOS = 0
  DO JM=0,R%NSMAX
    IF (ISENDSET == D%NPROCM(JM)) THEN
      INENTR = (D%NLATLE(MYSETW)-D%NLATLS(MYSETW)+1)*(R%NTMAX-JM+2)
      IF (IPOS + INENTR < NCOMBFLEN) THEN
        DO JGL=D%NLATLS(MYSETW),D%NLATLE(MYSETW)
          JGLLOC = JGL - D%NLATLS(MYSETW) + 1
          DO JN=1,R%NTMAX-JM+2
            IPOS = IPOS + 1
            ZCOMBUF(IPOS) = PNM(D%NPMG(JM)+JN,JGLLOC)
          ENDDO
        ENDDO
      ELSE
        DO JGL=D%NLATLS(MYSETW),D%NLATLE(MYSETW)
          JGLLOC = JGL - D%NLATLS(MYSETW) + 1
          DO JN=1,R%NTMAX-JM+2
            IPOS = IPOS + 1
            ZCOMBUF(IPOS) = PNM(D%NPMG(JM)+JN,JGLLOC)
            IF (IPOS == NCOMBFLEN) THEN
              CALL MPL_SEND(zcombuf(1:ipos),KDEST=NPRCIDS(ISEND), &
               & KTAG=ITAG,CDSTRING='SUTRLE:')
              IPOS = 0
              ITAG = ITAG + 1
              LLEXACT = (JGL == D%NLATLE(MYSETW) .AND. JN == R%NTMAX-JM+2)
              IF (.NOT.LLEXACT) LLADMSG = .TRUE.
            ENDIF
          ENDDO
        ENDDO
      ENDIF
    ENDIF
  ENDDO

!*   send message (if not empty or if message has been split)

  IF (IPOS > 0 .OR. LLADMSG) THEN
    CALL MPL_SEND(ZCOMBUF(1:IPOS),KDEST=NPRCIDS(ISEND), &
     & KTAG=ITAG,CDSTRING='SUTRLE:')
  ENDIF
  CALL GSTATS(801,1)

  ILREC = 0
  ITAG = MTAGLETR
  IF (D%NUMP > 0.AND. D%NLATLE(IRECSET) >= D%NLATLS(IRECSET)) THEN

!*   receive message (if not empty)

    CALL GSTATS(801,0)
    CALL MPL_RECV(ZCOMBUF(1:NCOMBFLEN),KSOURCE=NPRCIDS(IRECV), &
     & KTAG=ITAG,KOUNT=ILREC,CDSTRING='SUTRLE:')

!*   copy data from buffer to f%rpnm

    IPOS = 0
    DO JMLOC=1,D%NUMP
      JM = D%MYMS(JMLOC)
      INENTR = (D%NLATLE(IRECSET)-D%NLATLS(IRECSET)+1)*(R%NTMAX-JM+2)
      IOFFT = D%NPMT(JM) 
      IF (IPOS + INENTR < NCOMBFLEN) THEN
        DO JGL=D%NLATLS(IRECSET),D%NLATLE(IRECSET)
          DO JN=1,R%NTMAX-JM+2
            IPOS = IPOS + 1
            ZPNM(JGL,IOFFT+JN) = ZCOMBUF(IPOS)
          ENDDO
        ENDDO
      ELSE
        DO JGL=D%NLATLS(IRECSET),D%NLATLE(IRECSET)
          DO JN=1,R%NTMAX-JM+2
            IPOS = IPOS + 1
            ZPNM(JGL,IOFFT+JN) = ZCOMBUF(IPOS)
            IF (IPOS == NCOMBFLEN) THEN
              ITAG = ITAG + 1
              CALL MPL_RECV(ZCOMBUF(1:NCOMBFLEN), &
               & KSOURCE=NPRCIDS(IRECV),KTAG=ITAG, &
               & KOUNT=ILREC,CDSTRING='SUTRLE:')
              IPOS = 0
            ENDIF
          ENDDO
        ENDDO
      ENDIF
    ENDDO
    CALL GSTATS(801,1)

!*    check received message length

    IF (ILREC /= IPOS) THEN
      WRITE(NOUT,*)' SUTRLE: ILREC,IPOS,NCOMBLEN ',ILREC,IPOS,NCOMBFLEN
      CALL ABORT_TRANS(' SUTRLE:RECEIVED MESSAGE LENGTH DOES NOT MATCH')
    ENDIF
  ENDIF

! Perform barrier synchronisation to guarantee all processors have
! completed communication for this jroc loop iteration

  CALL MPL_BARRIER(CDSTRING='SUTRLE:')

ENDDO

!*    copy data from pnm to rpnm

CALL GSTATS(1803,0)
!cjfe OMP not efficient in that case
!cjfe!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jmloc,im,iofft,ioffg,jgl,iglloc,jn)
DO JMLOC=1,D%NUMP
  IM = D%MYMS(JMLOC)
  IOFFT = D%NPMT(IM)
  IOFFG = D%NPMG(IM)
  DO JGL=D%NLATLS(MYSETW),D%NLATLE(MYSETW)
    IGLLOC = JGL-D%NLATLS(MYSETW)+1
    DO JN=1,R%NTMAX-IM+2
      ZPNM(JGL,IOFFT+JN) = PNM(IOFFG+JN,IGLLOC)
    ENDDO
  ENDDO
ENDDO
!cjfe!$OMP END PARALLEL DO
CALL GSTATS(1803,1)

DEALLOCATE (ZCOMBUF)

END SUBROUTINE SUTRLE
END MODULE SUTRLE_MOD
