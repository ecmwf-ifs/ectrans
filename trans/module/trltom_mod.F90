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
!     ------------------------------------------------------------------

#include "tsmbkind.h"
USE MPL_MODULE
USE YOMGSTATS, ONLY : LSYNCSTATS

USE TPM_DISTR
USE TPM_GEN

USE SET2PE_MOD
USE MYSENDSET_MOD
USE MYRECVSET_MOD
USE ABORT_TRANS_MOD

IMPLICIT NONE


INTEGER_M,INTENT(IN)  :: KFIELD
REAL_B   ,INTENT(INOUT)  :: PFBUF(:)
REAL_B   ,INTENT(INOUT)  :: PFBUF_IN(:)
! === END OF INTERFACE BLOCK ===
INTEGER_M :: INFS(NPRTRW)

INTEGER_M :: ISENT    (NPRTRW-1)
INTEGER_M :: IRCVD    (NPRTRW-1)
INTEGER_M :: ISENDTOT (NPRTRW-1)
INTEGER_M :: IRECVTOT (NPRTRW-1)
INTEGER_M :: ISENDREQ (NPRTRW-1)
INTEGER_M :: IRECVREQ (NPRTRW-1)
INTEGER_M :: IJPOS(NPROC), IJSTOP(NPROC),  IJSEND(NPROC)
INTEGER_M :: ILENS(NPRTRW),IOFFS(NPRTRW),ILENR(NPRTRW),IOFFR(NPRTRW)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IBUFLEN, ICOMBFLENP, IERR, IFIRS1, IFIRSK1,&
             &IFIRSK2, IFIRST, IK1, IK2, ILAST, ILAST1, &
             &ILEN, ILENOL, ILREC, IP, IP1, IPOS, IPOSB, &
             &IREC, IRECD, IRECSET, IRECV, IRECVSET, ISEND,ISTP, &
             &ISENDSET, ISTA, ISTAR, ITAG, ITREC, J, J1, &
             &JHASE, JK, JROCA,IJR,IJS,INUMSENT

!     LOCAL LOGICAL SCALARS
LOGICAL :: LLDONE, LLEXIST
!     INTEGER FUNCTIONS


!     ------------------------------------------------------------------


!     TRANSPOSITION
!     -------------

! Set maximum transfer length

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
  IF (.NOT.LSYNCSTATS) CALL GSTATS(806,0)
  CALL MPL_ALLTOALLV(PSENDBUF=PFBUF_IN,KSENDCOUNTS=ILENS,&
   & PRECVBUF=PFBUF,KRECVCOUNTS=ILENR,KSENDDISPL=IOFFS,KRECVDISPL=IOFFR,&
   & KCOMM=MPL_ALL_MS_COMM,CDSTRING='TRLTOM:')
  
  IF (.NOT.LSYNCSTATS) CALL GSTATS(806,1)

ELSE
! Copy local part of PFBUF_IN to PFBUF

  ILEN = D%NLTSGTB(MYSETW)*KFIELD
  ISTA = D%NSTAGT1B(MYSETW)*KFIELD+1
  IF (.NOT.LSYNCSTATS) CALL GSTATS(1607,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(J)
  DO J=ISTA,ISTA+ILEN-1
    PFBUF(J) = PFBUF_IN(J)
  ENDDO
!$OMP END PARALLEL DO
  IF (.NOT.LSYNCSTATS) CALL GSTATS(1607,1)

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
           & KCOUNT=ILREC,CDSTRING='TRLTOM: LIMP WAITS ' )
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


!     ------------------------------------------------------------------

END SUBROUTINE TRLTOM
END MODULE TRLTOM_MOD
