! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


SUBROUTINE ETRANS_RELEASE(KRESOL)

!**** *ETRANS_RELEASE* - release a spectral resolution

!     Purpose.
!     --------
!      Release all arrays related to a given resolution tag

!**   Interface.
!     ----------
!     CALL ETRANS_RELEASE

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

USE PARKIND1  ,ONLY : JPIM

!ifndef INTERFACE

USE EDEALLOC_RESOL_MOD   ,ONLY : EDEALLOC_RESOL
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KRESOL

!endif INTERFACE

!     ------------------------------------------------------------------

CALL EDEALLOC_RESOL(KRESOL)

!     ------------------------------------------------------------------

END SUBROUTINE ETRANS_RELEASE
