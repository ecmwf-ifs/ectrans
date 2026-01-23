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

MODULE FTDIR_MOD
  USE BUFFERED_ALLOCATOR_MOD ,ONLY : ALLOCATION_RESERVATION_HANDLE
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: FTDIR, FTDIR_HANDLE, PREPARE_FTDIR

  TYPE FTDIR_HANDLE
    TYPE(ALLOCATION_RESERVATION_HANDLE) :: HREEL_COMPLEX
  END TYPE
CONTAINS

  FUNCTION PREPARE_FTDIR(ALLOCATOR,KF_FS) RESULT(HFTDIR)
    USE PARKIND_ECTRANS,        ONLY: JPIM, JPRBT, JPIB
    USE TPM_DISTR,              ONLY: D
    USE BUFFERED_ALLOCATOR_MOD, ONLY: BUFFERED_ALLOCATOR, RESERVE
    USE ISO_C_BINDING,          ONLY: C_SIZEOF

    IMPLICIT NONE

    TYPE(BUFFERED_ALLOCATOR), INTENT(INOUT) :: ALLOCATOR
    INTEGER(KIND=JPIM), INTENT(IN) :: KF_FS
    TYPE(FTDIR_HANDLE) :: HFTDIR

    REAL(KIND=JPRBT) :: DUMMY

#ifndef IN_PLACE_FFT
    HFTDIR%HREEL_COMPLEX = RESERVE(ALLOCATOR, 1_JPIB*KF_FS*D%NLENGTF*C_SIZEOF(DUMMY), "HFTDIR%HREEL_COMPLEX")
#endif
  END FUNCTION PREPARE_FTDIR

  SUBROUTINE FTDIR(ALLOCATOR,HFTDIR,PREEL_REAL,PREEL_COMPLEX,KFIELD)
    !**** *FTDIR - Direct Fourier transform

    !     Purpose. Routine for Grid-point to Fourier transform
    !     --------

    !**   Interface.
    !     ----------
    !        CALL FTDIR(..)

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

    USE TPM_GEN,                ONLY: LSYNC_TRANS, NCUR_RESOL
    USE PARKIND_ECTRANS,        ONLY: JPIM, JPRBT, JPIB
    USE TPM_DISTR,              ONLY: MYSETW, NPROC, D
    USE TPM_GEOMETRY,           ONLY: G
    USE BUFFERED_ALLOCATOR_MOD, ONLY: BUFFERED_ALLOCATOR, ASSIGN_PTR, GET_ALLOCATION
    USE TPM_HICFFT,             ONLY: EXECUTE_DIR_FFT
    USE MPL_MODULE,             ONLY: MPL_BARRIER,MPL_ALL_MS_COMM
    USE TPM_STATS,              ONLY: GSTATS => GSTATS_NVTX
    USE ISO_C_BINDING,          ONLY: C_SIZEOF

    IMPLICIT NONE

    INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELD
    REAL(KIND=JPRBT), INTENT(INOUT), POINTER :: PREEL_REAL(:)
    REAL(KIND=JPRBT), INTENT(OUT), POINTER :: PREEL_COMPLEX(:)
    TYPE(BUFFERED_ALLOCATOR), INTENT(IN) :: ALLOCATOR
    TYPE(FTDIR_HANDLE) :: HFTDIR

    INTEGER(KIND=JPIM) :: KGL

#ifdef IN_PLACE_FFT
    PREEL_COMPLEX => PREEL_REAL
#else
    CALL ASSIGN_PTR(PREEL_COMPLEX, GET_ALLOCATION(ALLOCATOR, HFTDIR%HREEL_COMPLEX),&
      & 1_JPIB, 1_JPIB*KFIELD*D%NLENGTF*C_SIZEOF(PREEL_COMPLEX(1)))
#endif

    ASSOCIATE(D_NDGL_FS=>D%NDGL_FS, D_NSTAGT0B=>D%NSTAGT0B, D_NSTAGTF=>D%NSTAGTF, &
            & D_NPTRLS=>D%NPTRLS, D_NPNTGTB0=>D%NPNTGTB0, D_NPROCM=>D%NPROCM, &
            & G_NMEN=>G%NMEN, G_NLOEN=>G%NLOEN)

#ifdef ACCGPU
    !$ACC DATA PRESENT(PREEL_REAL, PREEL_COMPLEX, &
    !$ACC&             D_NSTAGTF,D_NSTAGT0B,D_NPTRLS,D_NPROCM,D_NPNTGTB0,G_NMEN,G_NLOEN)
#endif
#ifdef OMPGPU
    !$OMP TARGET DATA MAP(PRESENT,ALLOC:PREEL_COMPLEX, &
    !$OMP&                D_NSTAGTF,D_NSTAGT0B,D_NPTRLS,D_NPROCM,D_NPNTGTB0,G_NMEN,G_NLOEN)
#endif

    IF (LSYNC_TRANS) THEN
      CALL GSTATS(430,0)
      CALL MPL_BARRIER(MPL_ALL_MS_COMM,CDSTRING='')
      CALL GSTATS(430,1)
    ENDIF
    CALL GSTATS(413,0)
    CALL EXECUTE_DIR_FFT(PREEL_REAL,PREEL_COMPLEX,NCUR_RESOL,KFIELD, &
        & LOENS=G_NLOEN(D_NPTRLS(MYSETW):D_NPTRLS(MYSETW)+D_NDGL_FS-1), &
        & OFFSETS=D_NSTAGTF(1:D_NDGL_FS+1),ALLOC=ALLOCATOR%PTR)

    IF (LSYNC_TRANS) THEN
      CALL GSTATS(433,0)
      CALL MPL_BARRIER(MPL_ALL_MS_COMM,CDSTRING='')
      CALL GSTATS(433,1)
    ENDIF
    CALL GSTATS(413,1)

#ifdef ACCGPU
    !$ACC END DATA
#endif
#ifdef OMPGPU
    !$OMP END TARGET DATA
#endif

    NULLIFY(PREEL_REAL)
    END ASSOCIATE

    !     ------------------------------------------------------------------
  END SUBROUTINE FTDIR
END MODULE FTDIR_MOD
