MODULE TRMTOL_MOD

CONTAINS
SUBROUTINE TRMTOL(PFBUF_IN,PFBUF,KFIELD)

!**** *trmtol * - transposition in Fourier space

!     Purpose.
!     --------
!              Transpose Fourier buffer data from partitioning
!              over wave numbers to partitioning over latitudes.
!              It is called between direct FFT and direct Legendre 
!              transform.
!              This routine is the inverse of TRLTOM.


!**   Interface.
!     ----------
!        *call* *trmtol(...)*

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
!        Modified : 97-06-17 G. Mozdzynski - control MPI mailbox use
!                                            (NCOMBFLEN) for nphase.eq.1
!        Modified : 99-05-28  D.Salmond - Optimise copies. 
!        Modified : 00-02-02  M.Hamrud  - Remove NPHASE
!        D.Salmond : 01-11-23 LIMP_NOOLAP Option for non-overlapping message
!                             passing and buffer packing
!     ------------------------------------------------------------------


#include "tsmbkind.h"
USE MPL_MODULE

USE TPM_DISTR
USE TPM_GEN

USE SET2PE_MOD
USE MYSENDSET_MOD
USE MYRECVSET_MOD
USE ABORT_TRANS_MOD

IMPLICIT NONE


INTEGER_M,INTENT(IN)    :: KFIELD
REAL_B   ,INTENT(INOUT) :: PFBUF(:)
REAL_B   ,INTENT(INOUT) :: PFBUF_IN(:)
! === END OF INTERFACE BLOCK ===

INTEGER_M :: ISENT    (NPRTRW-1)
INTEGER_M :: IRCVD    (NPRTRW-1)
INTEGER_M :: ISENDTOT (NPRTRW-1)
INTEGER_M :: IRECVTOT (NPRTRW-1)
INTEGER_M :: ISENDREQ (NPRTRW-1)
INTEGER_M :: IRECVREQ (NPRTRW-1)
INTEGER_M :: ILENS(NPRTRW),IOFFS(NPRTRW),ILENR(NPRTRW),IOFFR(NPRTRW)

!     LOCAL INTEGER SCALARS
INTEGER_M :: ICOMBFLENP, IERR,&
             &ILEN,  ILREC, &
             &IRECD, IRECV, IRECVSET, ISEND, ISTP,&
             &ISENDSET, ISTA, ITAG, ITREC, J, INUMSENT

!     LOCAL LOGICAL SCALARS
LOGICAL :: LLDONE
!     INTEGER FUNCTIONS


!     ------------------------------------------------------------------


!      transposition for nphase  = 1
!      -----------------------------


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
  ISENDTOT(J) = D%NLTSFTB(ISENDSET)*KFIELD
  IRECVTOT(J) = D%NLTSGTB(IRECVSET)*KFIELD
ENDDO

LLDONE = .FALSE.
ITAG = MTAGML



! Loop over the number of processors we need to communicate with

IF(LIMP_NOOLAP)THEN

!  Send/Recv loop.........................................................

  DO J=1,NPRTRW
    ILENS(J) = D%NLTSFTB(J)*KFIELD
    IOFFS(J) = D%NSTAGT0B(J)*KFIELD
    ILENR(J) = D%NLTSGTB(J)*KFIELD
    IOFFR(J) = D%NSTAGT0B(D%MSTABF(J))*KFIELD
  ENDDO

  CALL GSTATS(706,0)
  CALL MPL_ALLTOALLV(PSENDBUF=PFBUF_IN,KSENDCOUNTS=ILENS,&
   & PRECVBUF=PFBUF,KRECVCOUNTS=ILENR,KSENDDISPL=IOFFS,KRECVDISPL=IOFFR,&
   & KCOMM=MPL_ALL_MS_COMM,CDSTRING='TRMTOL:')
  CALL GSTATS(706,1)

ELSE

! For MYPROC copy straight accross

  ILEN = D%NLTSGTB(MYSETW)*KFIELD
  ISTA = D%NSTAGT0B(MYSETW)*KFIELD+1
  CALL GSTATS(1708,0)
!$OMP PARALLEL DO PRIVATE(J)
  DO J=ISTA,ISTA+ILEN-1
    PFBUF(J) = PFBUF_IN(J)
  ENDDO
!$OMP END PARALLEL DO
  CALL GSTATS(1708,1)

  DO WHILE( .NOT.LLDONE )
    
    LLDONE = .TRUE.
      
  
! Send loop.............................................................
  
!  For immediate send/recv, post receives first
  
    IF(LIMP) THEN
      DO J=1,NPRTRW-1
        IF( IRECVTOT(J)-IRCVD(J) > 0 )THEN
          IRECVSET = MYRECVSET(NPRTRW,MYSETW,J)
          CALL SET2PE(IRECV,0,0,IRECVSET,MYSETV)
          ILEN = MIN(ICOMBFLENP,IRECVTOT(J)-IRCVD(J))
          ISTA = D%NSTAGT0B(D%MSTABF(IRECVSET))*KFIELD+IRCVD(J)+1
          ISTP = ISTA+ILEN-1
          CALL MPL_RECV(PFBUF(ISTA:ISTP),KSOURCE=NPRCIDS(IRECV),KTAG=ITAG, &
           & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IRECVREQ(J), &
           & CDSTRING='TRMTOL:' )
        ENDIF
      ENDDO
    ENDIF
  
    INUMSENT = 0
  
    DO J=1,NPRTRW-1
  
! Check if there is no outstanding receive, and we have data to send
  
      IF(ISENDTOT(J)-ISENT(J) > 0 )THEN
        LLDONE = .FALSE.
        ISENDSET = MYSENDSET(NPRTRW,MYSETW,J)
        CALL SET2PE(ISEND,0,0,ISENDSET,MYSETV)
        ILEN = MIN(ICOMBFLENP,ISENDTOT(J)-ISENT(J))
        ISTA = D%NSTAGT0B(ISENDSET)*KFIELD+ISENT(J)+1
        ISTP = ISTA+ILEN-1
        IF(LIMP) THEN
          INUMSENT = INUMSENT+1
          CALL MPL_SEND(PFBUF_IN(ISTA:ISTP),KDEST=NPRCIDS(ISEND),&
           & KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(INUMSENT), &
           & KTAG=ITAG,CDSTRING='TRMTOL:')
        ELSE
          CALL MPL_SEND(PFBUF_IN(ISTA:ISTP),KDEST=NPRCIDS(ISEND),&
           &KTAG=ITAG,CDSTRING='TRMTOL:')
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
          ISTA = D%NSTAGT0B(D%MSTABF(IRECVSET))*KFIELD+IRCVD(J)+1
          ISTP = ISTA+ILEN-1
          CALL MPL_RECV(PFBUF(ISTA:ISTP),KSOURCE=NPRCIDS(IRECV),&
           &KTAG=ITAG,KOUNT=ILREC,CDSTRING='TRMTOL:')
        ELSE
!       For LIMP=true, find message length
          CALL MPL_WAIT(PFBUF,KREQUEST=IRECVREQ(J), &
           & KCOUNT=ILREC,CDSTRING='TRMTOL: LIMP WAITS ' )
        ENDIF
        IF( ILREC /= ILEN )THEN
          CALL ABORT_TRANS('TRMTOL:RECEIVED MESSAGE LENGTH ERROR')
        ENDIF
        IRCVD(J) = IRCVD(J)+ILEN
      ENDIF
    ENDDO
  
!       For LIMP=true, wait for sends to complete
    IF(LIMP) THEN
      IF(INUMSENT > 0) THEN
        CALL MPL_WAIT(PFBUF,KREQUEST=ISENDREQ(1:INUMSENT), &
         & CDSTRING='TRLTOM: ERROR IN MPL_WAIT FOR SENDS')
      ENDIF
    ENDIF
  
  ENDDO

! Perform barrier synchronisation to guarantee all processors have
! completed communication

  IF( NPRTRW > 1 )THEN
    CALL MPL_BARRIER(CDSTRING='TRMTOL:')
  ENDIF
ENDIF


!     ------------------------------------------------------------------

END SUBROUTINE TRMTOL
END MODULE TRMTOL_MOD
