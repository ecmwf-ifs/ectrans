! (C) Copyright 2025- Meteo-France.
! (C) Copyright 2025- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

SUBROUTINE ECTRANS_VERSION(CD_VERSION_STRING)
! ** PURPOSE
!    Return the version string of ecTrans
!
! ** DUMMY ARGUMENTS
!    CD_VERSION_STRING: version string
!
! ** AUTHOR
!    18 March 2025, S. Hatfield
!
! I. Dummy arguments declaration
USE ECTRANS_VERSION_MOD, ONLY: ECTRANS_VERSION_STR
USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR
IMPLICIT NONE
CHARACTER(KIND=C_CHAR), DIMENSION(14), INTENT(OUT) :: CD_VERSION_STRING
!
! II. Local variables
INTEGER :: JI
CHARACTER(LEN=SIZE(CD_VERSION_STRING)) :: C_VERSION_STRING
!
! III. Get version
C_VERSION_STRING = ECTRANS_VERSION_STR()
DO JI=1, SIZE(CD_VERSION_STRING)
  CD_VERSION_STRING(JI)=C_VERSION_STRING(JI:JI)
ENDDO
!
END SUBROUTINE ECTRANS_VERSION
