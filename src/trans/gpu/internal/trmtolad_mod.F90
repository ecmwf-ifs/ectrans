! (C) Copyright 1995- ECMWF.
! (C) Copyright 1995- Meteo-France.
! (C) Copyright 2022- NVIDIA.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TRMTOL_MOD
  USE BUFFERED_ALLOCATOR_MOD, ONLY: ALLOCATION_RESERVATION_HANDLE
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: TRMTOL, PREPARE_TRMTOL, TRMTOL_HANDLE

  TYPE TRMTOL_HANDLE
    TYPE(ALLOCATION_RESERVATION_HANDLE) :: HPFBUF
  END TYPE
CONTAINS
  FUNCTION PREPARE_TRMTOL(ALLOCATOR, KF_LEG) RESULT(HTRMTOL)
    USE PARKIND_ECTRANS,        ONLY: JPIM, JPRBT, JPIB
    USE TPM_DISTR,              ONLY: D
    USE BUFFERED_ALLOCATOR_MOD, ONLY: BUFFERED_ALLOCATOR, RESERVE
    USE ISO_C_BINDING,          ONLY: C_SIZE_T, C_SIZEOF

    IMPLICIT NONE

    TYPE(BUFFERED_ALLOCATOR), INTENT(INOUT) :: ALLOCATOR
    INTEGER(KIND=JPIM), INTENT(IN) :: KF_LEG
    TYPE(TRMTOL_HANDLE) :: HTRMTOL

    REAL(KIND=JPRBT) :: DUMMY

    HTRMTOL%HPFBUF = RESERVE(ALLOCATOR, 2_JPIB*D%NLENGT0B*KF_LEG*C_SIZEOF(DUMMY), "HTRMTOL%HPFBUF")
  END FUNCTION

  SUBROUTINE TRMTOL(ALLOCATOR,HTRMTOL,PFBUF_IN,PFBUF,KF_LEG)
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
    !                             KF_LEG - Number of fields communicated

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

    USE PARKIND_ECTRANS,        ONLY: JPIM, JPRBT, JPIB
    USE YOMHOOK,                ONLY: LHOOK, DR_HOOK, JPHOOK
    USE MPL_MODULE,             ONLY: MPL_ALLTOALLV, MPL_BARRIER, MPL_ALL_MS_COMM, MPL_MYRANK
    USE TPM_DISTR,              ONLY: D, NPRTRW, NPROC, MYPROC, MYSETW
    USE TPM_GEN,                ONLY: LSYNC_TRANS, NERR, LMPOFF
#if ECTRANS_HAVE_MPI
    USE MPI_F08,                ONLY: MPI_COMM, MPI_REAL4, MPI_REAL8
    ! Missing: MPI_ALLTOALLV on purpose due to cray-mpi bug (see https://github.com/ecmwf-ifs/ectrans/pull/157)
#endif
    USE BUFFERED_ALLOCATOR_MOD, ONLY: BUFFERED_ALLOCATOR, ASSIGN_PTR, GET_ALLOCATION
    USE TPM_STATS,              ONLY: GSTATS => GSTATS_NVTX
    USE ISO_C_BINDING,          ONLY: C_SIZE_T, C_SIZEOF
    USE ABORT_TRANS_MOD,        ONLY: ABORT_TRANS

    IMPLICIT NONE

    INTEGER(KIND=JPIM) ,INTENT(IN)  :: KF_LEG
    REAL(KIND=JPRBT), INTENT(OUT), POINTER  :: PFBUF(:)
    REAL(KIND=JPRBT), INTENT(IN) :: PFBUF_IN(:)

    INTEGER(KIND=JPIM) :: ILENS(NPRTRW),IOFFS(NPRTRW),ILENR(NPRTRW),IOFFR(NPRTRW)
    INTEGER(KIND=JPIM) :: J, FROM_SEND, TO_SEND, FROM_RECV, TO_RECV, IRANK
    INTEGER(KIND=JPIB) :: JPOS, ISTA, IEND, ILEN
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    INTEGER(KIND=JPIM) :: IERROR

    TYPE(BUFFERED_ALLOCATOR), INTENT(IN) :: ALLOCATOR
    TYPE(TRMTOL_HANDLE), INTENT(IN) :: HTRMTOL

#if ECTRANS_HAVE_MPI
    TYPE(MPI_COMM) :: LOCAL_COMM
#endif

#ifdef PARKINDTRANS_SINGLE
#define TRMTOL_DTYPE MPI_REAL4
#else
#define TRMTOL_DTYPE MPI_REAL8
#endif

#if ECTRANS_HAVE_MPI
    IF(.NOT. LMPOFF) THEN
      LOCAL_COMM%MPI_VAL = MPL_ALL_MS_COMM
    ENDIF
#endif

    IF (LHOOK) CALL DR_HOOK('TRMTOL',0,ZHOOK_HANDLE)

    CALL ASSIGN_PTR(PFBUF, GET_ALLOCATION(ALLOCATOR, HTRMTOL%HPFBUF),&
        & 1_JPIB, 2_JPIB*D%NLENGT0B*KF_LEG*C_SIZEOF(PFBUF(1)))

    IF(NPROC > 1) THEN
      DO J=1,NPRTRW
        ILENS(J) = D%NLTSFTB(J)*2*KF_LEG
        IOFFS(J) = D%NSTAGT1B(J)*2*KF_LEG
        ILENR(J) = D%NLTSGTB(J)*2*KF_LEG
        IOFFR(J) = D%NSTAGT0B(J)*2*KF_LEG
      ENDDO

      CALL GSTATS(807,0)

      ! copy to self workaround
      IRANK = MPL_MYRANK(MPL_ALL_MS_COMM)
      IF (ILENS(IRANK) .ne. ILENR(IRANK)) THEN
          WRITE(NERR,*) "ERROR", ILENS(IRANK), ILENR(IRANK)
          CALL ABORT_TRANS("TRMTOL: ILENS(IRANK) /= ILENR(IRANK)")
      ENDIF
      IF (ILENS(IRANK) > 0) THEN
          FROM_SEND = IOFFS(IRANK) + 1
          TO_SEND = FROM_SEND + ILENS(IRANK) - 1
          FROM_RECV = IOFFR(IRANK) + 1
          TO_RECV = FROM_RECV + ILENR(IRANK) - 1
#ifdef OMPGPU
#endif
#ifdef ACCGPU
#ifdef __HIP_PLATFORM_AMD__
          ! Workaround for AMD GPUs - ASYNC execution of this kernel gives numerical errors
          !$ACC KERNELS DEFAULT(NONE) PRESENT(PFBUF,PFBUF_IN) COPYIN(FROM_RECV,TO_RECV,FROM_SEND,TO_SEND)
#else
          !$ACC KERNELS ASYNC(1) DEFAULT(NONE) PRESENT(PFBUF,PFBUF_IN) COPYIN(FROM_RECV,TO_RECV,FROM_SEND,TO_SEND)
#endif
#endif
          PFBUF(FROM_RECV:TO_RECV) = PFBUF_IN(FROM_SEND:TO_SEND)
#ifdef OMPGPU
#endif
#ifdef ACCGPU
          !$ACC END KERNELS
#endif
          ILENS(IRANK) = 0
          ILENR(IRANK) = 0
      ENDIF

      IF (LSYNC_TRANS) THEN
        CALL GSTATS(440,0)
        CALL MPL_BARRIER(MPL_ALL_MS_COMM,CDSTRING='')
        CALL GSTATS(440,1)
      ENDIF
      CALL GSTATS(421,0)
#ifdef USE_GPU_AWARE_MPI
#ifdef OMPGPU
#endif
#ifdef ACCGPU
      !$ACC HOST_DATA USE_DEVICE(PFBUF_IN, PFBUF)
#endif
#else
    !! this is safe-but-slow fallback for running without GPU-aware MPI
    !$ACC UPDATE HOST(PFBUF_IN,PFBUF)
#endif

#if ECTRANS_HAVE_MPI
      CALL MPI_ALLTOALLV(PFBUF_IN,ILENS,IOFFS,TRMTOL_DTYPE,&
       & PFBUF,ILENR,IOFFR,TRMTOL_DTYPE,&
       & LOCAL_COMM,IERROR)
#else
      CALL ABORT_TRANS("Should not be here: MPI is disabled")
#endif

#ifdef USE_GPU_AWARE_MPI
#ifdef OMPGPU
#endif
#ifdef ACCGPU
      !$ACC END HOST_DATA
#endif
#else
    !! this is safe-but-slow fallback for running without GPU-aware MPI
    !$ACC UPDATE DEVICE(PFBUF)
#endif
      IF (LSYNC_TRANS) THEN
        CALL GSTATS(441,0)
        CALL MPL_BARRIER(MPL_ALL_MS_COMM,CDSTRING='')
        CALL GSTATS(441,1)
      ENDIF
      CALL GSTATS(421,1)

#ifdef ACCGPU
#ifndef __HIP_PLATFORM_AMD__
      ! Workaround for AMD GPUs - ASYNC execution of this kernel gives numerical errors
      !$ACC WAIT(1)
#endif
#endif
      CALL GSTATS(807,1)
    ELSE
      ILEN = 2_JPIB*D%NLTSGTB(MYSETW)*KF_LEG
      ISTA = 2_JPIB*D%NSTAGT0B(MYSETW)*KF_LEG+1
      IEND = ISTA+ILEN-1
      CALL GSTATS(1608,0)
#ifdef OMPGPU
#endif
#ifdef ACCGPU
      !$ACC PARALLEL LOOP DEFAULT(NONE) PRESENT(PFBUF,PFBUF_IN) FIRSTPRIVATE(ISTA,IEND)
#endif
      DO JPOS=ISTA,IEND
        PFBUF(JPOS) = PFBUF_IN(JPOS)
      ENDDO
      CALL GSTATS(1608,1)
    ENDIF

    IF (LHOOK) CALL DR_HOOK('TRMTOL',1,ZHOOK_HANDLE)

    !     ------------------------------------------------------------------
  END SUBROUTINE TRMTOL
END MODULE TRMTOL_MOD
