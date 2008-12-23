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
!        G.Mozdzynski : 08-01-01 Cleanup
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

INTEGER(KIND=JPIM) :: ILENS(NPRTRW),IOFFS(NPRTRW),ILENR(NPRTRW),IOFFR(NPRTRW)

INTEGER(KIND=JPIM) :: ITAG, J, ILEN, ISTA

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TRLTOM',0,ZHOOK_HANDLE)


ITAG = MTAGLM

DO J=1,NPRTRW
  ILENS(J) = D%NLTSGTB(J)*KFIELD
  IOFFS(J) = D%NSTAGT1B(D%MSTABF(J))*KFIELD
  ILENR(J) = D%NLTSFTB(J)*KFIELD
  IOFFR(J) = D%NSTAGT1B(J)*KFIELD
ENDDO

IF(NPROC > 1) THEN
  CALL GSTATS_BARRIER(763)
  CALL GSTATS(806,0)
  IF (LSYNC_TRANS) THEN
    CALL MPL_BARRIER(CDSTRING='TRLTOM:')
  ENDIF
  
  CALL MPL_ALLTOALLV(PSENDBUF=PFBUF_IN,KSENDCOUNTS=ILENS,&
   & PRECVBUF=PFBUF,KRECVCOUNTS=ILENR,KSENDDISPL=IOFFS,KRECVDISPL=IOFFR,&
   & KCOMM=MPL_ALL_MS_COMM,CDSTRING='TRLTOM:')
  CALL GSTATS(806,1)
  CALL GSTATS_BARRIER2(763)
ELSE
  ILEN = D%NLTSGTB(MYSETW)*KFIELD
  ISTA = D%NSTAGT1B(MYSETW)*KFIELD+1
  CALL GSTATS(1607,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(J)
  DO J=ISTA,ISTA+ILEN-1
    PFBUF(J) = PFBUF_IN(J)
  ENDDO
!$OMP END PARALLEL DO
  CALL GSTATS(1607,1)
ENDIF

IF (LHOOK) CALL DR_HOOK('TRLTOM',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------
END SUBROUTINE TRLTOM
END MODULE TRLTOM_MOD
