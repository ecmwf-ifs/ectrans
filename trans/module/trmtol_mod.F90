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
!     ------------------------------------------------------------------


#include "tsmbkind.h"
USE MPL_MODULE

USE TPM_DISTR

USE SET2PE_MOD


IMPLICIT NONE


INTEGER_M,INTENT(IN)    :: KFIELD
REAL_B   ,INTENT(INOUT) :: PFBUF(:)
REAL_B   ,INTENT(INOUT) :: PFBUF_IN(:)
! === END OF INTERFACE BLOCK ===

LOGICAL LLWAITRECV(NPRTRW)
INTEGER_M :: ISENT    (NPRTRW-1)
INTEGER_M :: IRCVD    (NPRTRW-1)
INTEGER_M :: ISENDTOT (NPRTRW-1)
INTEGER_M :: IRECVTOT (NPRTRW-1)

!     LOCAL INTEGER SCALARS
INTEGER_M :: ICOMBFLENP, IERR,&
             &ILEN,  ILREC, &
             &IRECD, IRECV, IRECVSET, ISEND, ISTP,&
             &ISENDSET, ISTA, ITAG, ITREC, J

!     LOCAL LOGICAL SCALARS
LOGICAL :: LLDONE, LLEXIST, LLWAIT, LLRECV
!     INTEGER FUNCTIONS
INTEGER_M :: MYSENDSET,MYRECVSET


!     ------------------------------------------------------------------


!      transposition for nphase  = 1
!      -----------------------------


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
  ISENDTOT(J)=D%NLTSFTB(ISENDSET)*KFIELD
  IRECVTOT(J)=D%NLTSGTB(IRECVSET)*KFIELD
ENDDO

LLDONE=.FALSE.
ITAG=MTAGML

! For MYPROC copy straight accross

ILEN=D%NLTSGTB(MYSETW)*KFIELD
ISTA=D%NSTAGT0B(MYSETW)*KFIELD+1

!$OMP PARALLEL DO PRIVATE(J)
DO J=ISTA,ISTA+ILEN-1
  PFBUF(J)=PFBUF_IN(J)
ENDDO
!$OMP END PARALLEL DO

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
      ISTA=D%NSTAGT0B(ISENDSET)*KFIELD+ISENT(J)+1
      ISTP=ISTA+ILEN-1
      CALL MPL_SEND(PFBUF_IN(ISTA:ISTP),KDEST=NPRCIDS(ISEND),&
       &KTAG=ITAG,CDSTRING='TRMTOL:')

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

    LLWAIT=.TRUE.
    CALL MPL_PROBE(KTAG=ITAG,LDWAIT=LLWAIT,LDFLAG=LLEXIST,&
     & CDSTRING='TRMTOL: PROBEING FOR ANY MESSAGE ')

  ENDIF

! A message should now be ready to be received

  DO J=1,NPRTRW-1
    IRECVSET=MYRECVSET(NPRTRW,MYSETW,J)
    IF( IRECVTOT(J)-IRCVD(J) > 0 )THEN
      CALL SET2PE(IRECV,0,0,IRECVSET,MYSETV)
      LLWAIT=.FALSE.
      CALL MPL_PROBE(KSOURCE=NPRCIDS(IRECV),KTAG=ITAG,LDWAIT=LLWAIT,&
       &LDFLAG=LLEXIST,CDSTRING='TRMTOL: PROBEING FOR MESSAGE ')
      IF( LLEXIST )THEN
        ILEN=MIN(ICOMBFLENP,IRECVTOT(J)-IRCVD(J))
        ISTA=D%NSTAGT0B(D%MSTABF(IRECVSET))*KFIELD+IRCVD(J)+1
        ISTP=ISTA+ILEN-1
        CALL MPL_RECV(PFBUF(ISTA:ISTP),KSOURCE=NPRCIDS(IRECV),&
         &KTAG=ITAG,KOUNT=ILREC,CDSTRING='TRMTOL:')
        IF( ILREC /= ILEN )THEN
          CALL ABOR1('TRMTOL:RECEIVED MESSAGE LENGTH ERROR')
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
  CALL MPL_BARRIER(CDSTRING='TRMTOL(END):')
ENDIF

RETURN
END SUBROUTINE TRMTOL
END MODULE TRMTOL_MOD
