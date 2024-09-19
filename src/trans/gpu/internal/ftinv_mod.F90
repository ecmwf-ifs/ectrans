! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! (C) Copyright 2022- NVIDIA.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FTINV_MOD
  USE BUFFERED_ALLOCATOR_MOD ,ONLY : BUFFERED_ALLOCATOR, ALLOCATION_RESERVATION_HANDLE
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: FTINV, FTINV_HANDLE, PREPARE_FTINV

  TYPE FTINV_HANDLE
    TYPE(ALLOCATION_RESERVATION_HANDLE) :: HREEL_REAL
  END TYPE
CONTAINS
  FUNCTION PREPARE_FTINV(ALLOCATOR,KF_FS) RESULT(HFTINV)
    USE PARKIND_ECTRANS,        ONLY: JPIM, JPRBT
    USE TPM_DISTR, ONLY: D
    USE BUFFERED_ALLOCATOR_MOD, ONLY: BUFFERED_ALLOCATOR, RESERVE
    USE ISO_C_BINDING,          ONLY: C_SIZE_T

    IMPLICIT NONE

    TYPE(BUFFERED_ALLOCATOR), INTENT(INOUT) :: ALLOCATOR
    INTEGER(KIND=JPIM), INTENT(IN) :: KF_FS
    TYPE(FTINV_HANDLE) :: HFTINV

    REAL(KIND=JPRBT) :: DUMMY

#ifndef IN_PLACE_FFT
    HFTINV%HREEL_REAL = RESERVE(ALLOCATOR, int(D%NLENGTF*KF_FS*SIZEOF(DUMMY),kind=c_size_t))
#endif
  END FUNCTION

  SUBROUTINE FTINV(ALLOCATOR,HFTINV,PREEL_COMPLEX,PREEL_REAL,KFIELD)
    !**** *FTINV - Inverse Fourier transform

    !     Purpose. Routine for Fourier to Grid-point transform
    !     --------

    !**   Interface.
    !     ----------
    !        CALL FTINV(..)

    !        Explicit arguments :  PREEL   - Fourier/grid-point array
    !        --------------------  KFIELD   - number of fields

    !     Method.
    !     -------

    !     Externals.  FFT992 - FFT routine
    !     ----------
    !

    !     Author.
    !     -------
    !        Mats Hamrud *ECMWF*

    !     Modifications.
    !     --------------
    !        Original : 00-03-03
    !        G. Radnoti 01-04-24 2D model (NLOEN=1)
    !        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
    !        G. Mozdzynski (Oct 2014): support for FFTW transforms
    !        G. Mozdzynski (Jun 2015): Support alternative FFTs to FFTW
    !     ------------------------------------------------------------------

    USE TPM_GEN,                ONLY: LSYNC_TRANS
    USE PARKIND_ECTRANS,        ONLY: JPIM, JPRBT
    USE TPM_DISTR,              ONLY: MYSETW, D_NPTRLS, D_NDGL_FS, D_NSTAGTF, D
    USE TPM_GEOMETRY,           ONLY: G_NLOEN
    USE TPM_HICFFT,             ONLY: EXECUTE_INV_FFT
    USE MPL_MODULE,             ONLY: MPL_BARRIER,MPL_ALL_MS_COMM
    USE TPM_STATS,              ONLY: GSTATS => GSTATS_NVTX
    USE BUFFERED_ALLOCATOR_MOD, ONLY: ASSIGN_PTR, GET_ALLOCATION
    USE ISO_C_BINDING,          ONLY: C_SIZE_T

    IMPLICIT NONE

    INTEGER(KIND=JPIM),INTENT(IN) :: KFIELD
    REAL(KIND=JPRBT), INTENT(INOUT), POINTER :: PREEL_REAL(:)
    REAL(KIND=JPRBT), INTENT(OUT), POINTER :: PREEL_COMPLEX(:)
    TYPE(BUFFERED_ALLOCATOR), INTENT(IN) :: ALLOCATOR
    TYPE(FTINV_HANDLE), INTENT(IN) :: HFTINV

    INTEGER(KIND=JPIM) :: KGL

#ifdef IN_PLACE_FFT
    PREEL_REAL => PREEL_COMPLEX
#else
    CALL ASSIGN_PTR(PREEL_REAL, GET_ALLOCATION(ALLOCATOR, HFTINV%HREEL_REAL),&
      & 1_C_SIZE_T, int(KFIELD*D%NLENGTF*SIZEOF(PREEL_REAL(1)),kind=c_size_t))
#endif

#ifdef OMPGPU
#endif
#ifdef ACCGPU
    !$ACC DATA PRESENT(PREEL_REAL,PREEL_COMPLEX,D_NPTRLS,D_NDGL_FS,D_NSTAGTF,G_NLOEN)
#endif

    IF (LSYNC_TRANS) THEN
      CALL GSTATS(440,0)
      CALL MPL_BARRIER(MPL_ALL_MS_COMM,CDSTRING='')
      CALL GSTATS(440,1)
    ENDIF
    CALL GSTATS(423,0)
    CALL EXECUTE_INV_FFT(PREEL_COMPLEX(:),PREEL_REAL(:),KFIELD, &
        & LOENS=G_NLOEN(D_NPTRLS(MYSETW):D_NPTRLS(MYSETW)+D_NDGL_FS-1), &
        & OFFSETS=D_NSTAGTF(1:D_NDGL_FS),ALLOC=ALLOCATOR%PTR)

    IF (LSYNC_TRANS) THEN
      CALL GSTATS(443,0)
      CALL MPL_BARRIER(MPL_ALL_MS_COMM,CDSTRING='')
      CALL GSTATS(443,1)
    ENDIF
    CALL GSTATS(423,1)

#ifdef OMPGPU
#endif
#ifdef ACCGPU
    !$ACC END DATA
#endif

    NULLIFY(PREEL_COMPLEX)

    !     ------------------------------------------------------------------
  END SUBROUTINE FTINV
END MODULE FTINV_MOD
