! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE TRANS_RELEASE(KRESOL)

!**** *TRANS_RELEASE* - release a spectral resolution

!     Purpose.
!     --------
!      Release all arrays related to a given resolution tag

!**   Interface.
!     ----------
!     CALL TRANS_RELEASE

!     Explicit arguments : KRESOL : resolution tag
!     --------------------

!     Method.
!     -------

!     Externals.  None
!     ----------

!     Author.
!     -------
!        R. El Khatib *METEO-FRANCE*

!     Modifications.
!     --------------
!        Original : 09-Jul-2013

!     ------------------------------------------------------------------

USE EC_PARKIND, ONLY: JPIM

!ifndef INTERFACE

USE RESOLS_MOD, ONLY: Y_RESOLS
USE TPM_GEN, ONLY: LENABLED, NOUT

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KRESOL

!endif INTERFACE

!     ------------------------------------------------------------------

! Check this resol is less than maximum number of resols
IF (KRESOL < Y_RESOLS%SIZE) THEN
  ! Check this resol is actually initialised
  IF (LENABLED(KRESOL)) THEN
    CALL Y_RESOLS(KRESOL)%DESTROY(KRESOL)
  ELSE
    WRITE(NOUT, FMT='(''TRANS_RELEASE: Warning KRESOL = '',I3,'' is already disabled '')') KRESOL
  ENDIF
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE TRANS_RELEASE
