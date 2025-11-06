! (C) Copyright 2015- ECMWF.
! (C) Copyright 2015- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE VD2UV_MOD
CONTAINS
SUBROUTINE VD2UV(KM,KMLOC,KF_UV,KLEI2,PSPVOR,PSPDIV,PU,PV)

USE PARKIND_ECTRANS, ONLY: JPIM, JPRB
USE YOMHOOK,         ONLY: LHOOK, DR_HOOK, JPHOOK
USE ABORT_TRANS_MOD, ONLY: ABORT_TRANS


!**** *VD2UV* - U and V from Vor/div
!
!     Purpose.
!     --------
!       
!**   Interface.
!     ----------
!        *CALL* *VD2UV(...)

!        Explicit arguments :
!        --------------------
!          KM        - zonal wavenumber
!          KMLOC     - local zonal wavenumber
!          PSPVOR    - spectral vorticity
!          PSPDIV    - spectral divergence
!          PU(:,:)   - spectral U (out)
!          PV(:,:)   - spectral V (out)


!        Implicit arguments :  

!     Method.
!     -------

!     Externals.
!     ----------

!         PREPSNM - prepare REPSNM for wavenumber KM
!         PRFI1B  - prepares the spectral fields
!         VDTUV   - compute u and v from vorticity and divergence

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Temperton, 1991, MWR 119 p1303

!     Author.
!     -------
!        Mats Hamrud  *ECMWF*

!     Modifications.
!     --------------
!        Original : July 2015
!       
!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KM
INTEGER(KIND=JPIM), INTENT(IN) :: KMLOC
INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV
INTEGER(KIND=JPIM), INTENT(IN) :: KLEI2

REAL(KIND=JPRB)   , INTENT(IN)  :: PSPVOR(:,:)
REAL(KIND=JPRB)   , INTENT(IN)  :: PSPDIV(:,:)
REAL(KIND=JPRB)   , INTENT(OUT) :: PU(:,:)
REAL(KIND=JPRB)   , INTENT(OUT) :: PV(:,:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------

IF (LHOOK) CALL DR_HOOK('VD2UV_MOD',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       1.    PREPARE ZEPSNM.
!              ---------------

CALL ABORT_TRANS('VD2UV: Code path not (yet) supported in GPU version')

IF (LHOOK) CALL DR_HOOK('VD2UV_MOD',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE VD2UV
END MODULE VD2UV_MOD




