! (C) Copyright 2014- ECMWF.
! (C) Copyright 2014- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE CDMAP_MOD
CONTAINS
SUBROUTINE CDMAP(KM,KMLOC,KSL,KSLO,PEPSNM, KDIR, KDGNH, KDGNHD,&
&                KFIELDS, PCOEFA, PCOEFS)

USE PARKIND_ECTRANS, ONLY: JPIM,  JPRB, JPRBT
USE MPL_MODULE,      ONLY: MPL_ABORT

!**** *CDMAP* - REMAP ROOTS
!
!     Purpose.
!     --------
! remap from one set of roots to another using Christoffel-Darboux formula, see Chien + Alpert, 1997.

!**   Interface.
!     ----------
!        *CALL* *CDMAP(...)

!        Explicit arguments :
!        --------------------
!          KM        - zonal wavenumber
!          KMLOC     - local zonal wavenumber
!
!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!      Chien + Alpert, 1997.

!     Author.
!     -------
!        Nils Wedi  *ECMWF*

!     Modifications.
!     --------------
!        Original : 14-05-14 
!     ------------------------------------------------------------------

IMPLICIT NONE


INTEGER(KIND=JPIM), INTENT(IN) :: KM
INTEGER(KIND=JPIM), INTENT(IN) :: KMLOC
INTEGER(KIND=JPIM), INTENT(IN) :: KSL
INTEGER(KIND=JPIM), INTENT(IN) :: KSLO
REAL(KIND=JPRBT),   INTENT(IN) :: PEPSNM
INTEGER(KIND=JPIM), INTENT(IN) :: KDIR ! direction of map
INTEGER(KIND=JPIM), INTENT(IN) :: KDGNH
INTEGER(KIND=JPIM), INTENT(IN) :: KDGNHD
INTEGER(KIND=JPIM), INTENT(IN) :: KFIELDS
REAL(KIND=JPRBT),   INTENT(INOUT) :: PCOEFA(:,:)
REAL(KIND=JPRBT),   INTENT(INOUT) :: PCOEFS(:,:)

!     ------------------------------------------------------------------

CALL MPL_ABORT("CDMAP not yet supported in ecTrans GPU version")

!     ------------------------------------------------------------------

END SUBROUTINE CDMAP
END MODULE CDMAP_MOD
