! (C) Copyright 2026- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE RESOL_CHECK4PY_MOD
! Validation of the resolution tag KRESOL of the distributed ectrans4py routines, so that
! a tag which is not a resolution set up by SETUP_TRANS_4PY comes back to Python as an
! error code (RuntimeError) instead of ending the whole MPI job: ecTrans's SET_RESOL calls
! ABORT_TRANS for it. The check is local but gives the same result on every task, since
! all tasks set up the same resolutions.
!
! ecTrans offers no public inquiry of whether a resolution is set up, so RESOL_STATUS
! reads the state SET_RESOL checks, in the internal module TPM_GEN: NMAX_RESOL (the
! maximum number of resolutions, KMAX_RESOL of SETUP_TRANS0) and LENABLED (the
! resolutions set up and not released).
USE ISO_FORTRAN_ENV, ONLY: INT64
IMPLICIT NONE
PRIVATE
PUBLIC :: RESOL_STATUS, JPRESOL_UNKNOWN

! Return code for a KRESOL that is not a resolution set up by SETUP_TRANS_4PY
INTEGER(KIND=INT64), PARAMETER :: JPRESOL_UNKNOWN = 101_INT64

CONTAINS

FUNCTION RESOL_STATUS(KRESOL) RESULT(KSTATUS)
! 0 if KRESOL is a resolution set up by SETUP_TRANS_4PY, else JPRESOL_UNKNOWN (with a
! message on the error unit).
USE ISO_FORTRAN_ENV, ONLY: ERROR_UNIT
USE TPM_GEN, ONLY: NMAX_RESOL, LENABLED
INTEGER(KIND=INT64), INTENT(IN) :: KRESOL
INTEGER(KIND=INT64) :: KSTATUS
KSTATUS = 0_INT64
IF (KRESOL < 1_INT64 .OR. KRESOL > INT(NMAX_RESOL, INT64)) THEN
  KSTATUS = JPRESOL_UNKNOWN
ELSEIF (.NOT. ALLOCATED(LENABLED)) THEN
  KSTATUS = JPRESOL_UNKNOWN
ELSEIF (.NOT. LENABLED(KRESOL)) THEN
  KSTATUS = JPRESOL_UNKNOWN
ENDIF
IF (KSTATUS /= 0_INT64) THEN
  WRITE(ERROR_UNIT, '(A,I0,A,I0,A)') 'ectrans4py: KRESOL = ', KRESOL, &
   & ' is not a resolution set up by SETUP_TRANS_4PY (maximum number of resolutions: ', &
   & NMAX_RESOL, ')'
ENDIF
END FUNCTION RESOL_STATUS

END MODULE RESOL_CHECK4PY_MOD
