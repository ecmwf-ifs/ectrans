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

MODULE SET_RESOL_MOD
CONTAINS
SUBROUTINE SET_RESOL(KRESOL,LDSETUP)
USE PARKIND1,        ONLY: JPIM, JPRD, JPRB
USE TPM_GEN,         ONLY: NOUT, MSETUP0, NCUR_RESOL, NMAX_RESOL, LENABLED
USE TPM_DIM,         ONLY: R, DIM_RESOL
USE TPM_DISTR,       ONLY: D, DISTR_RESOL
USE TPM_GEOMETRY,    ONLY: G, GEOM_RESOL
USE TPM_FIELDS,      ONLY: F, FIELDS_RESOL
USE TPM_FIELDS_GPU,  ONLY: FG, FIELDS_GPU_RESOL
USE TPM_FLT,         ONLY: S, FLT_RESOL
USE TPM_CTL,         ONLY: C, CTL_RESOL
USE ABORT_TRANS_MOD, ONLY: ABORT_TRANS
USE RESOLS_MOD,      ONLY: Y_RESOLS
use BACKENDS_MOD,    ONLY: JP_BACKEND_GPU_SP, JP_BACKEND_GPU_DP, JP_UNINITIALISED

IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KRESOL
LOGICAL            ,OPTIONAL, INTENT(IN) :: LDSETUP

! Local variables
INTEGER(KIND=JPIM) :: IRESOL, I_BACKEND
LOGICAL :: LLSETUP

!     ------------------------------------------------------------------

IF(MSETUP0 == 0) CALL ABORT_TRANS('SET_RESOL:TRANS NOT SETUP')
LLSETUP = .FALSE.
IF(PRESENT(LDSETUP)) LLSETUP = LDSETUP
IRESOL = 1
IF(PRESENT(KRESOL)) THEN
  IRESOL = KRESOL
 IF(IRESOL < 1 .OR. IRESOL > NMAX_RESOL) THEN
   WRITE(NOUT,*)'SET_RESOL: UNKNOWN RESOLUTION ',IRESOL,NMAX_RESOL
   CALL ABORT_TRANS('SET_RESOL:IRESOL < 1 .OR. KRESOL > NMAX_RESOL')
  ENDIF
  IF(.NOT.LLSETUP) THEN
    IF(.NOT.LENABLED(IRESOL)) THEN
      WRITE(NOUT,*)'SET_RESOL: UNKNOWN RESOLUTION ',IRESOL,LENABLED
      CALL ABORT_TRANS('SET_RESOL:IRESOL NOT ENABLED')
     ENDIF
   ENDIF
ENDIF

! Verify backend compatibility
! etrans is not yet "backend-aware", but we don't currently prohibit the use of TRANS_INQ in
! combination with etrans, so for now we need to permit etrans to bypass this check.
! etrans is the only one that doesn't give Y_RESOLS(IRESOL)%BACKEND a value, so it will be
! JP_UNINITIALISED.
! In future it would be nice if we could enforce that ONLY ETRANS_INQ is to be used with etrans.
IF (Y_RESOLS(IRESOL)%BACKEND /= JP_UNINITIALISED) THEN
  I_BACKEND = MERGE(JP_BACKEND_GPU_DP, JP_BACKEND_GPU_SP, JPRB == JPRD)
  IF (Y_RESOLS(IRESOL)%BACKEND /= I_BACKEND) THEN
    WRITE(NOUT,*) 'SET_RESOL: RESOLUTION ', IRESOL, ' NOT SUPPORTED FOR BACKEND ', I_BACKEND
    CALL ABORT_TRANS('SET_RESOL:RESOLUTION NOT SUPPORTED FOR BACKEND')
  ENDIF
ENDIF

IF(IRESOL /= NCUR_RESOL) THEN
  NCUR_RESOL = IRESOL
ENDIF

! We always reassociate the pointers, even if IRESOL equalled NCUR_RESOL on entering this subroutine
! If a different backend called SET_RESOL with the same IRESOL previously, and that resol was
! deactivated, that IRESOL might be reused here, but the pointers will be pointing to the other
! backend's arrays, so they need to be reassociated
R => DIM_RESOL(NCUR_RESOL)
F => FIELDS_RESOL(NCUR_RESOL)
FG => FIELDS_GPU_RESOL(NCUR_RESOL)
G => GEOM_RESOL(NCUR_RESOL)
D => DISTR_RESOL(NCUR_RESOL)
S => FLT_RESOL(NCUR_RESOL)
C => CTL_RESOL(NCUR_RESOL)

END SUBROUTINE SET_RESOL
END MODULE SET_RESOL_MOD
