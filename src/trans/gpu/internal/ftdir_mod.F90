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
  IMPLICIT NONE

  TYPE FTDIR_HANDLE
  END TYPE
CONTAINS

  FUNCTION PREPARE_FTDIR() RESULT(HFTDIR)
    IMPLICIT NONE
    TYPE(FTDIR_HANDLE) :: HFTDIR
  END FUNCTION

  SUBROUTINE FTDIR(HFTDIR,PREEL_REAL,PREEL_COMPLEX,KFIELD)
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

    USE TPM_GEN         ,ONLY : LSYNC_TRANS
    USE PARKIND_ECTRANS ,ONLY : JPIM, JPRBT

    USE TPM_DISTR       ,ONLY : D, MYSETW, MYPROC, NPROC, D_NSTAGT0B, D_NSTAGTF,D_NPTRLS, D_NPNTGTB0, D_NPROCM
    USE TPM_GEOMETRY    ,ONLY : G, G_NMEN, G_NLOEN
    USE TPM_HICFFT      ,ONLY : EXECUTE_DIR_FFT
    USE MPL_MODULE      ,ONLY : MPL_BARRIER
    USE TPM_STATS, ONLY : GSTATS => GSTATS_NVTX

    IMPLICIT NONE

    INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELD
    REAL(KIND=JPRBT), INTENT(INOUT), POINTER :: PREEL_REAL(:)
    REAL(KIND=JPRBT), INTENT(OUT), POINTER :: PREEL_COMPLEX(:)
    TYPE(FTDIR_HANDLE) :: HFTDIR

    INTEGER(KIND=JPIM) :: KGL, D_NDGL_FS

    D_NDGL_FS = D%NDGL_FS

    PREEL_COMPLEX => PREEL_REAL

#ifdef ACCGPU
    !$ACC DATA PRESENT(PREEL_REAL, PREEL_COMPLEX, &
    !$ACC&             D_NSTAGTF,D_NSTAGT0B,D_NPTRLS,D_NPROCM,D_NPNTGTB0,G_NMEN,G_NLOEN)
#endif
#ifdef OMPGPU
    !$OMP TARGET DATA MAP(ALLOC:PREEL_COMPLEX)
#endif

    IF (LSYNC_TRANS) THEN
      CALL GSTATS(430,0)
      CALL MPL_BARRIER(CDSTRING='')
      CALL GSTATS(430,1)
    ENDIF
    CALL GSTATS(413,0)
    CALL EXECUTE_DIR_FFT(PREEL_REAL(:),PREEL_COMPLEX(:),KFIELD, &
        & LOENS=G%NLOEN(D%NPTRLS(MYSETW):D%NPTRLS(MYSETW)+D%NDGL_FS-1), &
        & OFFSETS=D%NSTAGTF(1:D%NDGL_FS+1))

    IF (LSYNC_TRANS) THEN
      CALL GSTATS(433,0)
      CALL MPL_BARRIER(CDSTRING='')
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

    !     ------------------------------------------------------------------
  END SUBROUTINE FTDIR
END MODULE FTDIR_MOD
