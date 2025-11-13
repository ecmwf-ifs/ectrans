! (C) Copyright 2025- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE SPECTRAL_SPACE_FUNC_MOD



IMPLICIT NONE

PRIVATE
PUBLIC SPECTRAL_SPACE_FUNC

INTERFACE
  SUBROUTINE SPECTRAL_SPACE_FUNC(KM, KSL, KDGL, KFIELDS, P1MU2, PFIELD, KPTRU, KFLDUV, &
    &                            KFLDSC, KFLDPTRUV)
  USE PARKIND1, ONLY: JPIM, JPRB
  INTEGER(KIND=JPIM), INTENT(IN) :: KM
  INTEGER(KIND=JPIM), INTENT(IN) :: KSL
  INTEGER(KIND=JPIM), INTENT(IN) :: KDGL
  INTEGER(KIND=JPIM), INTENT(IN) :: KFIELDS
  REAL(KIND=JPRB), INTENT(IN) :: P1MU2(KDGL)
  REAL(KIND=JPRB), INTENT(INOUT) :: PFIELD(2*KFIELDS,0:KDGL+1)
  INTEGER(KIND=JPIM), INTENT(IN) :: KPTRU
  INTEGER(KIND=JPIM), INTENT(IN) :: KFLDUV
  INTEGER(KIND=JPIM), INTENT(IN) :: KFLDSC
  INTEGER(KIND=JPIM), INTENT(IN) :: KFLDPTRUV(KFLDUV)
  END SUBROUTINE SPECTRAL_SPACE_FUNC
END INTERFACE

END MODULE SPECTRAL_SPACE_FUNC_MOD