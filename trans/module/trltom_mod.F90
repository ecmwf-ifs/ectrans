MODULE TRLTOM_MOD
CONTAINS
SUBROUTINE TRLTOM(PFBUF_IN,PFBUF,KFIELD)

!**** *trltom * - transposition in Fourierspace

!     Purpose.
!     --------
!              Transpose Fourier coefficients from partitioning
!              over latitudes to partitioning over wave numbers
!              This is done between inverse Legendre Transform
!              and inverse FFT.
!              This is the inverse routine of TRMTOL.

!**   Interface.
!     ----------
!        *call* *trltom(...)*

!        Explicit arguments : PFBUF  - Fourier coefficient buffer. It is
!        --------------------          used for both input and output.
!                                      KSTA 
!                             KFBUFL - Dimension of PFBUF
!                             KSTA   - Start address of data sent to/
!                                      received from other PE's
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
!     ------------------------------------------------------------------

#include "tsmbkind.h"

use tpm_distr

use set2pe_mod

IMPLICIT NONE


INTEGER_M,INTENT(IN)  :: KFIELD
REAL_B   ,INTENT(INOUT)  :: PFBUF(:)
REAL_B   ,INTENT(INOUT)  :: PFBUF_IN(:)
! === END OF INTERFACE BLOCK ===
INTEGER_M :: INFS(NPRTRW)

LOGICAL LLWAITRECV(NPRTRW)
INTEGER_M :: ISENT    (NPRTRW-1)
INTEGER_M :: IRCVD    (NPRTRW-1)
INTEGER_M :: ISENDTOT (NPRTRW-1)
INTEGER_M :: IRECVTOT (NPRTRW-1)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IBUFLEN, ICOMBFLENP, IERR, IFIRS1, IFIRSK1,&
             &IFIRSK2, IFIRST, IK1, IK2, ILAST, ILAST1, &
             &ILEN, ILENOL, ILREC, IP, IP1, IPOS, IPOSB, &
             &IREC, IRECD, IRECSET, IRECV, IRECVSET, ISEND, &
             &ISENDSET, ISTA, ISTAR, ITAG, ITREC, J, J1, &
             &JHASE, JK, JROCA

!     LOCAL LOGICAL SCALARS
LOGICAL :: LLDONE, LLEXIST, LLWAIT, LLRECV
!     INTEGER FUNCTIONS
INTEGER_M :: MYSENDSET,MYRECVSET


!     ------------------------------------------------------------------



!     transposition
!     -------------


! Set maximum transfer length

IF( NPRTRW > 1 )THEN
  ICOMBFLENP=NCOMBFLEN/(3*NPRTRW-4)
ENDIF

! Mark all processors as having no outstanding receives,
! nothing sent or received.

LLWAITRECV(:)=.FALSE.
DO J=1,NPRTRW-1
  ISENDSET=MYSENDSET(NPRTRW,MYSETW,J)
  IRECVSET=MYRECVSET(NPRTRW,MYSETW,J)
  ISENT(J)=0
  IRCVD(J)=0
  ISENDTOT(J)=D%NLTSGTB(ISENDSET)*KFIELD
  IRECVTOT(J)=D%NLTSFTB(IRECVSET)*KFIELD
ENDDO
LLDONE=.FALSE.
ITAG=MTAGLM

! Copy local part of PFBUF_IN to PFBUF

ILEN=D%NLTSGTB(MYSETW)*KFIELD
ISTA=D%NSTAGT1B(MYSETW)*KFIELD+1
PFBUF(ISTA:ISTA+ILEN-1)=PFBUF_IN(ISTA:ISTA+ILEN-1)


! Loop over the number of processors we need to communicate with

DO WHILE( .NOT.LLDONE )

  LLDONE=.TRUE.


! Send loop.............................................................


  DO J=1,NPRTRW-1

! Check if there is no outstanding receive, and we have data to send

    ISENDSET=MYSENDSET(NPRTRW,MYSETW,J)
    IF( .NOT.LLWAITRECV(ISENDSET).AND.ISENDTOT(J)-ISENT(J) > 0 )THEN

      LLDONE=.FALSE.
      CALL SET2PE(ISEND,0,0,ISENDSET,MYSETV)
      ILEN=MIN(ICOMBFLENP,ISENDTOT(J)-ISENT(J))
      ISTA=D%NSTAGT1B(D%MSTABF(ISENDSET))*KFIELD+ISENT(J)+1

      CALL MPE_SEND(PFBUF_IN(ISTA),ILEN,MREALT,NPRCIDS(ISEND),&
       &ITAG,0,0,0,IERR)
      IF( IERR < 0 )THEN
        CALL ABOR1('TRLTOM : ERROR IN MPE_SEND')
      ENDIF

      ISENT(J)=ISENT(J)+ILEN
      
    ENDIF

  ENDDO


!  Receive loop.........................................................


! Check if there is data to receive

  LLRECV=.FALSE.
  DO J=1,NPRTRW-1
    IF( IRECVTOT(J)-IRCVD(J) > 0 )THEN
      LLRECV=.TRUE.
      LLDONE=.FALSE.
    ENDIF
  ENDDO

  IF( LLRECV )THEN

! Probe and wait for a message from anybody

    IRECV=-1
    LLWAIT=.TRUE.
    CALL MPE_PROBE(IRECV,ITAG,0,LLWAIT,LLEXIST,IERR)
    IF( IERR < 0 )THEN
      CALL ABOR1('TRLTOM: ERROR IN MPE_PROBE WAIT=T')
    ENDIF

  ENDIF

! A message should now be ready to be received

  DO J=1,NPRTRW-1
    IRECVSET=MYRECVSET(NPRTRW,MYSETW,J)
    IF( IRECVTOT(J)-IRCVD(J) > 0 )THEN
      CALL SET2PE(IRECV,0,0,IRECVSET,MYSETV)
      LLWAIT=.FALSE.
      CALL MPE_PROBE(NPRCIDS(IRECV),ITAG,0,LLWAIT,LLEXIST,IERR)
      IF( IERR < 0 )THEN
        CALL ABOR1('TRLTOM: ERROR IN MPE_PROBE WAIT=F')
      ENDIF
      IF( LLEXIST )THEN
        ILEN=MIN(ICOMBFLENP,IRECVTOT(J)-IRCVD(J))
        ISTA=D%NSTAGT1B(IRECVSET)*KFIELD+IRCVD(J)+1
        CALL MPE_RECV(PFBUF(ISTA),ILEN,MREALT,NPRCIDS(IRECV),&
         &ITAG,0,0,0,ILREC,IRECD,ITREC,IERR)
        IF( IERR < 0 )THEN
          CALL ABOR1('TRLTOM : ERROR IN MPE_RECV')
        ENDIF
        IF( ILREC /= ILEN )THEN
          CALL ABOR1('TRLTOM:RECEIVED MESSAGE LENGTH ERROR')
        ENDIF
        IRCVD(J)=IRCVD(J)+ILEN
        LLWAITRECV(IRECVSET)=.FALSE.
      ELSE
! Expecting data but none found, set outstanding receive flag
        LLWAITRECV(IRECVSET)=.TRUE.
      ENDIF
    ELSE
! No data to receive
      LLWAITRECV(IRECVSET)=.FALSE.
    ENDIF
  ENDDO

ENDDO

! Perform barrier synchronisation to guarantee all processors have
! completed communication

IF( NPROC > 1 )THEN
  CALL MPE_BARRIER(IERR)
  IF( IERR /= 0 )THEN
    CALL ABOR1('TRLTOM: ERROR IN MPE_BARRIER')
  ENDIF
ENDIF

END SUBROUTINE TRLTOM
END MODULE TRLTOM_MOD
