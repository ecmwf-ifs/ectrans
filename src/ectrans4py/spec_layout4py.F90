! (C) Copyright 2026- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE SPEC_LAYOUT4PY_MOD
! Repacking of spectral arrays between the layout of the Python callers and that of ecTrans,
! for every binding with spectral fields.
!
! A C-contiguous NumPy array (nfld, kspec2) is, in Fortran storage order, (KSPEC2, KFLD), while
! INV_TRANS, DIR_TRANS and SPECNORM take the fields first, (KFLD, KSPEC2), so the bindings
! transpose their spectral arguments and results. There are two implementations:
!  * by default, the array-section loops
!      DO JL = 1, KFLD
!        POUT(JL, :) = PIN(:, JL)
!      ENDDO
!    in which the accesses to one of the two arrays are strided by KSPEC2;
!  * with the environment variable ECTRANS4PY_OMP set to a value that does not start with '0'
!    (e.g. 1), loops blocked over the spectral coefficients, JPBLK at a time, so that the strided
!    accesses stay within a tile in the cache, and with the blocks shared out among the OpenMP
!    threads of the task (OMP_NUM_THREADS). Whether they are faster than the array-section
!    loops depends on the compiler, the machine and the number of threads (with one compiler,
!    they were measured slower on a single thread), so they are not the default.
! Both are pure permutations of the values: their results are bit-identical, whatever the number
! of threads.
!
! This is the only source of ectrans4py compiled with OpenMP (see CMakeLists.txt), so that the
! code generated for the other sources is the same with and without OpenMP.
USE PARKIND1, ONLY: JPIM, JPRB
IMPLICIT NONE
PRIVATE
PUBLIC :: SPEC_TO_FLD, FLD_TO_SPEC, SPEC_LAYOUT_INQ

! Number of spectral coefficients in a block of the blocked loops
INTEGER(KIND=JPIM), PARAMETER :: JPBLK = 64
! Whether ECTRANS4PY_OMP selects the blocked loops, read from the environment on first use (threads
! racing the first use read and store the same value)
LOGICAL, SAVE :: LLBLOCKED = .FALSE.
LOGICAL, SAVE :: LLBLOCKED_READ = .FALSE.

CONTAINS

LOGICAL FUNCTION LL_BLOCKED()
! Whether the environment variable ECTRANS4PY_OMP is set to a value that does not start with '0'
CHARACTER(LEN=1) :: CLVAL
INTEGER :: ILEN, ISTAT
IF (.NOT. LLBLOCKED_READ) THEN
  CLVAL = ' '
  CALL GET_ENVIRONMENT_VARIABLE('ECTRANS4PY_OMP', CLVAL, ILEN, ISTAT)
  ! (status -1: set to a value longer than CLVAL, which holds its first character)
  LLBLOCKED = (ISTAT == 0 .OR. ISTAT == -1) .AND. ILEN >= 1 .AND. CLVAL /= '0'
  LLBLOCKED_READ = .TRUE.
ENDIF
LL_BLOCKED = LLBLOCKED
END FUNCTION LL_BLOCKED

SUBROUTINE SPEC_TO_FLD(KSPEC2, KFLD, PIN, POUT)
! POUT(KFLD, KSPEC2) = TRANSPOSE(PIN(KSPEC2, KFLD)): from the layout of the Python callers to
! that of ecTrans
INTEGER(KIND=JPIM), INTENT(IN)  :: KSPEC2, KFLD
REAL(KIND=JPRB),    INTENT(IN)  :: PIN(KSPEC2, KFLD)
REAL(KIND=JPRB),    INTENT(OUT) :: POUT(KFLD, KSPEC2)
INTEGER(KIND=JPIM) :: JB, JS, JL
IF (LL_BLOCKED()) THEN
  ! each thread writes whole columns of POUT; the strided reads of a block stay in the cache
!$OMP PARALLEL DO PRIVATE(JB, JS, JL) SCHEDULE(STATIC)
  DO JB = 1, KSPEC2, JPBLK
    DO JS = JB, MIN(JB + JPBLK - 1, KSPEC2)
      DO JL = 1, KFLD
        POUT(JL, JS) = PIN(JS, JL)
      ENDDO
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
ELSE
  DO JL = 1, KFLD
    POUT(JL, :) = PIN(:, JL)
  ENDDO
ENDIF
END SUBROUTINE SPEC_TO_FLD

SUBROUTINE FLD_TO_SPEC(KSPEC2, KFLD, PIN, POUT)
! POUT(KSPEC2, KFLD) = TRANSPOSE(PIN(KFLD, KSPEC2)): from the layout of ecTrans to that of the
! Python callers
INTEGER(KIND=JPIM), INTENT(IN)  :: KSPEC2, KFLD
REAL(KIND=JPRB),    INTENT(IN)  :: PIN(KFLD, KSPEC2)
REAL(KIND=JPRB),    INTENT(OUT) :: POUT(KSPEC2, KFLD)
INTEGER(KIND=JPIM) :: JB, JS, JL
IF (LL_BLOCKED()) THEN
  ! each thread writes a block of every column of POUT, reading whole columns of PIN
!$OMP PARALLEL DO PRIVATE(JB, JS, JL) SCHEDULE(STATIC)
  DO JB = 1, KSPEC2, JPBLK
    DO JL = 1, KFLD
      DO JS = JB, MIN(JB + JPBLK - 1, KSPEC2)
        POUT(JS, JL) = PIN(JL, JS)
      ENDDO
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
ELSE
  DO JL = 1, KFLD
    POUT(:, JL) = PIN(JL, :)
  ENDDO
ENDIF
END SUBROUTINE FLD_TO_SPEC

SUBROUTINE SPEC_LAYOUT_INQ(KBLK, LDBLOCKED, KTHREADS, LDOPENMP, LDOPENMP_BUILD)
! The block size KBLK of the blocked loops; whether ECTRANS4PY_OMP selects them (LDBLOCKED); the
! number of OpenMP threads they use (KTHREADS, 1 without OpenMP); whether this module is
! compiled with OpenMP (LDOPENMP), and whether the build asked for it (LDOPENMP_BUILD: the
! feature OMP of ecTrans). The two must agree, or the blocked loops would run on a single thread
! without anything else showing it.
!$ USE OMP_LIB, ONLY: OMP_GET_MAX_THREADS
INTEGER(KIND=JPIM), INTENT(OUT) :: KBLK, KTHREADS
LOGICAL,            INTENT(OUT) :: LDBLOCKED, LDOPENMP, LDOPENMP_BUILD
KBLK = JPBLK
LDBLOCKED = LL_BLOCKED()
KTHREADS = 1
LDOPENMP = .FALSE.
!$ KTHREADS = OMP_GET_MAX_THREADS()
!$ LDOPENMP = .TRUE.
#ifdef ECTRANS4PY_WITH_OPENMP
LDOPENMP_BUILD = .TRUE.
#else
LDOPENMP_BUILD = .FALSE.
#endif
END SUBROUTINE SPEC_LAYOUT_INQ

END MODULE SPEC_LAYOUT4PY_MOD


! Bindings of the repacking and its inquiry, so that they can be tested and timed from Python

SUBROUTINE SPEC_TO_FLD4PY(KRETURNCODE, KSPEC2, KFLD, PSPEC, PFLD)
! PFLD(KFLD, KSPEC2), the transpose of PSPEC(KSPEC2, KFLD), as the bindings with spectral fields
! repack their inputs (SPEC_TO_FLD)
USE ISO_FORTRAN_ENV, ONLY: INT64
USE PARKIND1, ONLY: JPIM, JPRB
USE SPEC_LAYOUT4PY_MOD, ONLY: SPEC_TO_FLD
IMPLICIT NONE
INTEGER(KIND=INT64), INTENT(OUT) :: KRETURNCODE
INTEGER(KIND=INT64), INTENT(IN)  :: KSPEC2, KFLD
REAL(KIND=JPRB), DIMENSION(KSPEC2, KFLD), INTENT(IN)  :: PSPEC
REAL(KIND=JPRB), DIMENSION(KFLD, KSPEC2), INTENT(OUT) :: PFLD
KRETURNCODE = 0
CALL SPEC_TO_FLD(INT(KSPEC2,JPIM), INT(KFLD,JPIM), PSPEC, PFLD)
END SUBROUTINE SPEC_TO_FLD4PY


SUBROUTINE FLD_TO_SPEC4PY(KRETURNCODE, KSPEC2, KFLD, PFLD, PSPEC)
! PSPEC(KSPEC2, KFLD), the transpose of PFLD(KFLD, KSPEC2), as the bindings with spectral fields
! repack their outputs (FLD_TO_SPEC)
USE ISO_FORTRAN_ENV, ONLY: INT64
USE PARKIND1, ONLY: JPIM, JPRB
USE SPEC_LAYOUT4PY_MOD, ONLY: FLD_TO_SPEC
IMPLICIT NONE
INTEGER(KIND=INT64), INTENT(OUT) :: KRETURNCODE
INTEGER(KIND=INT64), INTENT(IN)  :: KSPEC2, KFLD
REAL(KIND=JPRB), DIMENSION(KFLD, KSPEC2), INTENT(IN)  :: PFLD
REAL(KIND=JPRB), DIMENSION(KSPEC2, KFLD), INTENT(OUT) :: PSPEC
KRETURNCODE = 0
CALL FLD_TO_SPEC(INT(KSPEC2,JPIM), INT(KFLD,JPIM), PFLD, PSPEC)
END SUBROUTINE FLD_TO_SPEC4PY


SUBROUTINE SPEC_LAYOUT_INQ4PY(KRETURNCODE, KBLK, KBLOCKED, KTHREADS, KOPENMP, KOPENMP_BUILD)
! SPEC_LAYOUT_INQ, with the logicals as integers (1: true)
USE ISO_FORTRAN_ENV, ONLY: INT64
USE PARKIND1, ONLY: JPIM
USE SPEC_LAYOUT4PY_MOD, ONLY: SPEC_LAYOUT_INQ
IMPLICIT NONE
INTEGER(KIND=INT64), INTENT(OUT) :: KRETURNCODE
INTEGER(KIND=INT64), INTENT(OUT) :: KBLK, KBLOCKED, KTHREADS, KOPENMP, KOPENMP_BUILD
INTEGER(KIND=JPIM) :: IBLK, ITHREADS
LOGICAL :: LLBLOCKED, LLOPENMP, LLOPENMP_BUILD
KRETURNCODE = 0
CALL SPEC_LAYOUT_INQ(IBLK, LLBLOCKED, ITHREADS, LLOPENMP, LLOPENMP_BUILD)
KBLK = IBLK
KTHREADS = ITHREADS
KBLOCKED = MERGE(1_INT64, 0_INT64, LLBLOCKED)
KOPENMP = MERGE(1_INT64, 0_INT64, LLOPENMP)
KOPENMP_BUILD = MERGE(1_INT64, 0_INT64, LLOPENMP_BUILD)
END SUBROUTINE SPEC_LAYOUT_INQ4PY
