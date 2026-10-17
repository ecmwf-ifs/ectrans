! (C) Copyright 2026- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE RESOL_CHECK4PY_MOD
! Validation of the arguments of the distributed ectrans4py routines against the setup of
! ecTrans, so that a bad one comes back to Python as an error code (RuntimeError) instead of
! ending the whole MPI job or reading memory that does not exist:
!  * a resolution tag KRESOL that is not a resolution set up by SETUP_TRANS_4PY (ecTrans's
!    SET_RESOL calls ABORT_TRANS for it): error code 101;
!  * a grid-point-only resolution (SETUP_TRANS_4PY LDGRIDONLY) for a routine that needs its
!    spectral or Fourier space, which a grid-point-only resolution does not have (the
!    transforms, the spectral scatter and gather, the norms -- GPNORM_TRANS too, whose
!    TRGTOL transposition uses the Fourier-space layout -- and the inquiries of spectral
!    sizes): error code 102;
!  * the region counts of the partition tables of TRANS_INQ_GRID4PY (REGION_COUNTS).
! The checks are local but give the same result on every task, since all tasks set up the
! same resolutions.
!
! ecTrans offers no public inquiry of this state, so this module reads internal modules of
! ecTrans (the distributed bindings read ecTrans internals nowhere else):
!  * TPM_GEN: NMAX_RESOL (the maximum number of resolutions, KMAX_RESOL of SETUP_TRANS0)
!    and LENABLED (the resolutions set up and not released), which SET_RESOL checks;
!  * TPM_DISTR: DISTR_RESOL(KRESOL)%LGRIDONLY, the LDGRIDONLY of SETUP_TRANS;
!  * EQ_REGIONS_MOD: N_REGIONS_NS and N_REGIONS_EW, the region counts of SETUP_TRANS0.
USE ISO_FORTRAN_ENV, ONLY: INT64
IMPLICIT NONE
PRIVATE
PUBLIC :: RESOL_STATUS, JPRESOL_UNKNOWN, JPRESOL_GRIDONLY, REGION_COUNTS

! Return code for a KRESOL that is not a resolution set up by SETUP_TRANS_4PY
INTEGER(KIND=INT64), PARAMETER :: JPRESOL_UNKNOWN = 101_INT64
! Return code for a grid-point-only KRESOL, for a routine that needs spectral space
INTEGER(KIND=INT64), PARAMETER :: JPRESOL_GRIDONLY = 102_INT64

CONTAINS

FUNCTION RESOL_STATUS(KRESOL, LDGRIDPOINT_ONLY) RESULT(KSTATUS)
! 0 if the routine can work on resolution KRESOL, else JPRESOL_UNKNOWN or JPRESOL_GRIDONLY
! (with a message on the error unit). LDGRIDPOINT_ONLY (default .FALSE.): the routine uses
! grid-point space only, and so can work on a grid-point-only resolution too.
USE ISO_FORTRAN_ENV, ONLY: ERROR_UNIT
USE TPM_GEN, ONLY: NMAX_RESOL, LENABLED
USE TPM_DISTR, ONLY: DISTR_RESOL
INTEGER(KIND=INT64), INTENT(IN) :: KRESOL
LOGICAL, OPTIONAL, INTENT(IN) :: LDGRIDPOINT_ONLY
INTEGER(KIND=INT64) :: KSTATUS
LOGICAL :: LLGRIDPOINT_ONLY
LLGRIDPOINT_ONLY = .FALSE.
IF (PRESENT(LDGRIDPOINT_ONLY)) LLGRIDPOINT_ONLY = LDGRIDPOINT_ONLY
KSTATUS = 0_INT64
IF (KRESOL < 1_INT64 .OR. KRESOL > INT(NMAX_RESOL, INT64)) THEN
  KSTATUS = JPRESOL_UNKNOWN
ELSEIF (.NOT. ALLOCATED(LENABLED) .OR. .NOT. ALLOCATED(DISTR_RESOL)) THEN
  KSTATUS = JPRESOL_UNKNOWN
ELSEIF (.NOT. LENABLED(KRESOL)) THEN
  KSTATUS = JPRESOL_UNKNOWN
ELSEIF (DISTR_RESOL(KRESOL)%LGRIDONLY .AND. .NOT. LLGRIDPOINT_ONLY) THEN
  KSTATUS = JPRESOL_GRIDONLY
ENDIF
IF (KSTATUS == JPRESOL_UNKNOWN) THEN
  WRITE(ERROR_UNIT, '(A,I0,A,I0,A)') 'ectrans4py: KRESOL = ', KRESOL, &
   & ' is not a resolution set up by SETUP_TRANS_4PY (maximum number of resolutions: ', &
   & NMAX_RESOL, ')'
ELSEIF (KSTATUS == JPRESOL_GRIDONLY) THEN
  WRITE(ERROR_UNIT, '(A,I0,A)') 'ectrans4py: KRESOL = ', KRESOL, &
   & ' is a grid-point-only resolution (SETUP_TRANS_4PY LDGRIDONLY), without spectral space'
ENDIF
END FUNCTION RESOL_STATUS

SUBROUTINE REGION_COUNTS(KREGIONS_NS, KREGIONS_EW)
! The numbers of N-S bands and of E-W regions of the grid-point partition, as SETUP_TRANS0
! returns them (K_REGIONS_NS, K_REGIONS_EW), to check the sizes of the partition tables a
! caller passes to TRANS_INQ_GRID4PY. ecTrans offers no public inquiry of them, so they are
! read from the internal module EQ_REGIONS_MOD (N_REGIONS_NS, N_REGIONS_EW).
USE EQ_REGIONS_MOD, ONLY: N_REGIONS_NS, N_REGIONS_EW
INTEGER(KIND=INT64), INTENT(OUT) :: KREGIONS_NS, KREGIONS_EW
KREGIONS_NS = INT(N_REGIONS_NS, INT64)
KREGIONS_EW = INT(N_REGIONS_EW, INT64)
END SUBROUTINE REGION_COUNTS

END MODULE RESOL_CHECK4PY_MOD
