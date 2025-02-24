! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE EASRE1AD_MOD
CONTAINS
SUBROUTINE EASRE1AD(KM,KMLOC,KF_OUT_LT,PIA)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!USE TPM_TRANS
USE EASRE1BAD_MOD   ,ONLY : EASRE1BAD

!**** *EASRE1AD* - Recombine antisymmetric and symmetric parts - adjoint

!     Purpose.
!     --------
!        To recombine the antisymmetric and symmetric parts of the
!        Fourier arrays and update the correct parts of the state
!        variables.

!**   Interface.
!     ----------
!       *CALL* *EASRE1AD(...)

!        Explicit arguments :
!        --------------------
!                              KM - zonal wavenumber
!                              KMLOC - local zonal wavenumber
!                              PAOA1 - antisymmetric part of Fourier
!                              fields for zonal wavenumber KM (basic
!                              variables and N-S derivatives)
!                              PSOA1 - symmetric part of Fourier
!                              fields for zonal wavenumber KM (basic
!                              variables and N-S derivatives)

!        Implicit arguments :  None
!        --------------------

!     Method.
!     -------

!     Externals.   EASRE1BAD - basic recombination routine
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From ASRE1AD in IFS CY22R1
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM) , INTENT(IN)  :: KM
INTEGER(KIND=JPIM) , INTENT(IN)  :: KMLOC
INTEGER(KIND=JPIM) , INTENT(IN)  :: KF_OUT_LT

REAL(KIND=JPRB)    , INTENT(OUT) :: PIA(:,:)

INTEGER(KIND=JPIM) :: IFLDS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EASRE1AD_MOD:EASRE1AD',0,ZHOOK_HANDLE)
IFLDS = KF_OUT_LT

CALL EASRE1BAD(IFLDS,KM,KMLOC,PIA)
IF (LHOOK) CALL DR_HOOK('EASRE1AD_MOD:EASRE1AD',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE EASRE1AD
END MODULE EASRE1AD_MOD
