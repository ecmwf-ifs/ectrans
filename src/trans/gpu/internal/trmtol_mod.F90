! (C) Copyright 1995- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TRMTOL_MOD

CONTAINS
SUBROUTINE TRMTOL_CUDAAWARE(PFBUF_IN,PFBUF,KFIELD)

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
!        G.Mozdzynski: 08-01-01 Cleanup
!        Y.Seity   : 07-08-31 add barrier synchronisation under LSYNC_TRANS
!     ------------------------------------------------------------------

USE PARKIND_ECTRANS  ,ONLY : JPIM     ,JPRBT
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE MPL_MODULE  ,ONLY : MPL_ALLTOALLV, MPL_BARRIER, MPL_ALL_MS_COMM, MPL_MYRANK
USE TPM_DISTR       ,ONLY : D, NPRTRW, NPROC, MYPROC, MYSETW
USE TPM_GEN         ,ONLY : LSYNC_TRANS
USE MPI
USE TPM_STATS, ONLY : GSTATS => GSTATS_NVTX

IMPLICIT NONE

INTEGER(KIND=JPIM) ,INTENT(IN)  :: KFIELD
REAL(KIND=JPRBT), INTENT(OUT), ALLOCATABLE  :: PFBUF(:)
REAL(KIND=JPRBT), INTENT(IN), ALLOCATABLE  :: PFBUF_IN(:)

INTEGER(KIND=JPIM) :: ILENS(NPRTRW),IOFFS(NPRTRW),ILENR(NPRTRW),IOFFR(NPRTRW)
INTEGER(KIND=JPIM) :: J, ILEN, ISTA, FROM_SEND, TO_SEND, FROM_RECV, TO_RECV, IRANK
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER(KIND=JPIM) :: IERROR

#ifdef PARKINDTRANS_SINGLE
#define TRMTOL_DTYPE MPI_REAL
#else
#define TRMTOL_DTYPE MPI_DOUBLE_PRECISION
#endif

IF (LHOOK) CALL DR_HOOK('TRMTOL_CUDAAWARE',0,ZHOOK_HANDLE)

ALLOCATE(PFBUF(D%NLENGT0B*KFIELD))
!$ACC ENTER DATA CREATE(PFBUF)

IF(NPROC > 1) THEN
  DO J=1,NPRTRW
    ILENS(J) = D%NLTSFTB(J)*KFIELD
    IOFFS(J) = D%NSTAGT1B(J)*KFIELD
    ILENR(J) = D%NLTSGTB(J)*KFIELD
    IOFFR(J) = D%NSTAGT0B(J)*KFIELD
  ENDDO

  CALL GSTATS(807,0)

  ! copy to self workaround
  IRANK = MPL_MYRANK(MPL_ALL_MS_COMM)
  IF (ILENS(IRANK) .ne. ILENR(IRANK)) THEN
      PRINT *, "ERROR", ILENS(IRANK), ILENR(IRANK)
      stop 1
  ENDIF
  IF (ILENS(IRANK) > 0) THEN
      FROM_SEND = IOFFS(IRANK) + 1
      TO_SEND = FROM_SEND + ILENS(IRANK) - 1
      FROM_RECV = IOFFR(IRANK) + 1
      TO_RECV = FROM_RECV + ILENR(IRANK) - 1
      !$ACC KERNELS ASYNC(1) DEFAULT(NONE) PRESENT(PFBUF,PFBUF_IN)
      PFBUF(FROM_RECV:TO_RECV) = PFBUF_IN(FROM_SEND:TO_SEND)
      !$ACC END KERNELS
      ILENS(IRANK) = 0
      ILENR(IRANK) = 0
  ENDIF

  IF (LSYNC_TRANS) THEN
    CALL GSTATS(440,0)
    CALL MPL_BARRIER(CDSTRING='')
    CALL GSTATS(440,1)
  ENDIF
  CALL GSTATS(421,0)
  !$ACC HOST_DATA USE_DEVICE(PFBUF_IN, PFBUF)
  CALL MPI_ALLTOALLV(PFBUF_IN,ILENS,IOFFS,TRMTOL_DTYPE,&
   & PFBUF,ILENR,IOFFR,TRMTOL_DTYPE,&
   & MPL_ALL_MS_COMM,IERROR)
  !$ACC END HOST_DATA
  IF (LSYNC_TRANS) THEN
    CALL GSTATS(441,0)
    CALL MPL_BARRIER(CDSTRING='')
    CALL GSTATS(441,1)
  ENDIF
  CALL GSTATS(421,1)

  !$ACC WAIT(1)
  CALL GSTATS(807,1)
ELSE
  ILEN = D%NLTSGTB(MYSETW)*KFIELD
  ISTA = D%NSTAGT0B(MYSETW)*KFIELD+1
  CALL GSTATS(1608,0)
  !$ACC PARALLEL LOOP DEFAULT(NONE) PRESENT(PFBUF,PFBUF_IN)
  DO J=ISTA,ISTA+ILEN-1
    PFBUF(J) = PFBUF_IN(J)
  ENDDO
  CALL GSTATS(1608,1)
ENDIF

IF (ALLOCATED(PFBUF_IN)) THEN
  !$ACC EXIT DATA DELETE(PFBUF_IN)
  DEALLOCATE(PFBUF_IN)
ENDIF

IF (LHOOK) CALL DR_HOOK('TRMTOL_CUDAAWARE',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------
END SUBROUTINE TRMTOL_CUDAAWARE

SUBROUTINE TRMTOL(PFBUF_IN,PFBUF,KFIELD)

!**** *TRMTOL * - transposition in Fourier space

!     Purpose.
!     --------
!              Transpose Fourier buffer data from partitioning
!              over wave numbers to partitioning over latitudes.
!              It is called between direct FFT and direct Legendre
!              transform.
!              This routine is the inverse of TRLTOM.


!**   Interface.
!     ----------
!        *CALL* *TRMTOL(...)*

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
!        G.Mozdzynski: 08-01-01 Cleanup
!        Y.Seity   : 07-08-31 Add barrier synchronisation under LSYNC_TRANS
!     ------------------------------------------------------------------

USE PARKIND_ECTRANS  ,ONLY : JPIM     ,JPRBT
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE MPL_MODULE  ,ONLY : MPL_ALLTOALLV, MPL_BARRIER, MPL_ALL_MS_COMM
USE TPM_DISTR       ,ONLY : D, NPRTRW, NPROC, MYPROC, MYSETW
USE TPM_GEN         ,ONLY : LSYNC_TRANS
USE MPI

IMPLICIT NONE

INTEGER(KIND=JPIM) ,INTENT(IN)  :: KFIELD
REAL(KIND=JPRBT), INTENT(OUT), ALLOCATABLE  :: PFBUF(:)
REAL(KIND=JPRBT), INTENT(IN), ALLOCATABLE  :: PFBUF_IN(:)

INTEGER(KIND=JPIM) :: ILENS(NPRTRW),IOFFS(NPRTRW),ILENR(NPRTRW),IOFFR(NPRTRW)
INTEGER(KIND=JPIM) :: J, ILEN, ISTA
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER(KIND=JPIM) :: IERROR

IF (LHOOK) CALL DR_HOOK('TRMTOL',0,ZHOOK_HANDLE)

ALLOCATE(PFBUF(D%NLENGT0B*KFIELD))
!$ACC ENTER DATA CREATE(PFBUF)

IF(NPROC > 1) THEN
  DO J=1,NPRTRW
    ILENS(J) = D%NLTSFTB(J)*KFIELD
    IOFFS(J) = D%NSTAGT1B(J)*KFIELD
    ILENR(J) = D%NLTSGTB(J)*KFIELD
    IOFFR(J) = D%NSTAGT0B(J)*KFIELD
  ENDDO

  CALL GSTATS(807,0)
  !$ACC UPDATE HOST(PFBUF_IN)

  CALL MPL_ALLTOALLV(PSENDBUF=PFBUF_IN,KSENDCOUNTS=ILENS,&
   & PRECVBUF=PFBUF,KRECVCOUNTS=ILENR,KSENDDISPL=IOFFS,KRECVDISPL=IOFFR,&
   & KCOMM=MPL_ALL_MS_COMM,CDSTRING='TRMTOL:')

  !$ACC UPDATE DEVICE(PFBUF)
  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
  CALL GSTATS(807,1)
ELSE
  ILEN = D%NLTSGTB(MYSETW)*KFIELD
  ISTA = D%NSTAGT0B(MYSETW)*KFIELD+1
  CALL GSTATS(1608,0)
  !$ACC PARALLEL LOOP DEFAULT(NONE) PRESENT(PFBUF,PFBUF_IN)
  DO J=ISTA,ISTA+ILEN-1
    PFBUF(J) = PFBUF_IN(J)
  ENDDO
  CALL GSTATS(1608,1)
ENDIF

!$ACC EXIT DATA DELETE(PFBUF_IN)
DEALLOCATE(PFBUF_IN)

IF (LHOOK) CALL DR_HOOK('TRMTOL',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------
END SUBROUTINE TRMTOL
END MODULE TRMTOL_MOD
