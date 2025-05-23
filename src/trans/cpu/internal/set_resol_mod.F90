! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
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
USE PARKIND1  ,ONLY : JPIM

USE TPM_GEN         ,ONLY : NOUT, MSETUP0, NCUR_RESOL, NMAX_RESOL,LENABLED
USE TPM_DIM         ,ONLY : R, DIM_RESOL
USE TPM_DISTR       ,ONLY : D, DISTR_RESOL
USE TPM_GEOMETRY    ,ONLY : G, GEOM_RESOL
USE TPM_FIELDS      ,ONLY : F, FIELDS_RESOL
USE TPM_FFTW        ,ONLY : TW, FFTW_RESOL
USE TPM_FLT         ,ONLY : S, FLT_RESOL
USE TPM_CTL         ,ONLY : C, CTL_RESOL
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
!

IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KRESOL
LOGICAL            ,OPTIONAL, INTENT(IN) :: LDSETUP

! Local variables
INTEGER(KIND=JPIM) :: IRESOL
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
IF(IRESOL /= NCUR_RESOL) THEN
  NCUR_RESOL = IRESOL
  R => DIM_RESOL(NCUR_RESOL)
  F => FIELDS_RESOL(NCUR_RESOL)
  G => GEOM_RESOL(NCUR_RESOL)
  D => DISTR_RESOL(NCUR_RESOL)
  TW => FFTW_RESOL(NCUR_RESOL)
  S => FLT_RESOL(NCUR_RESOL)
  C => CTL_RESOL(NCUR_RESOL)
ENDIF

END SUBROUTINE SET_RESOL
END MODULE SET_RESOL_MOD
