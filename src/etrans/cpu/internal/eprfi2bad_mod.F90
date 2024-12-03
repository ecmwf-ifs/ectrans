! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE EPRFI2BAD_MOD
CONTAINS
SUBROUTINE EPRFI2BAD(KFIELD,KM,KMLOC,PFFT)

!**** *EPRFI2BAD* - Prepare input work arrays for direct transform

!     Purpose.
!     --------
!        To extract the Fourier fields for a specific zonal wavenumber
!        and put them in an order suitable for the direct Legendre
!        tranforms, i.e. split into symmetric and anti-symmetric part.

!**   Interface.
!     ----------
!     *CALL* *EPRFI2BAD(..)

!        Explicit arguments :
!        -------------------   KFIELD - number of fields
!                              KM - zonal wavenumber
!                              KMLOC - local zonal wavenumber
!                              PAOA - antisymmetric part of Fourier
!                              fields for zonal wavenumber KM
!                              PSOA - symmetric part of Fourier
!                              fields for zonal wavenumber KM

!        Implicit arguments :  FOUBUF in TPM_TRANS
!        --------------------

!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 90-07-01
!        MPP Group: 95-10-01 Support for Distributed Memory version
!        Modified : 04/06/99 D.Salmond : change order of AIA and SIA
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DIM         ,ONLY : R
!USE TPMALD_DIM      ,ONLY : RALD
USE TPM_TRANS       ,ONLY : FOUBUF
!USE TPM_GEOMETRY
USE TPM_DISTR       ,ONLY : D
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELD,KM,KMLOC
REAL(KIND=JPRB)  , INTENT(IN)  :: PFFT(:,:)

INTEGER(KIND=JPIM) :: ISTAN, JF, JGL

INTEGER(KIND=JPIM) :: IJR,IJI
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.    EXTRACT SYM./ANTISYM. FIELDS FROM FOURIER ARRAY.
!              ------------------------------------------------
IF (LHOOK) CALL DR_HOOK('EPRFI2BAD_MOD:EPRFI2BAD',0,ZHOOK_HANDLE)
DO JGL=1,R%NDGL
  ISTAN = (D%NSTAGT1B(D%NPROCL(JGL) )+D%NPNTGTB1(KMLOC,JGL ))*2*KFIELD
  DO JF =1,KFIELD
    IJR = 2*(JF-1)+1
    IJI = IJR+1
    FOUBUF(ISTAN+IJR) = PFFT(JGL,IJR)
    FOUBUF(ISTAN+IJI) = PFFT(JGL,IJI)
  ENDDO
ENDDO
IF (LHOOK) CALL DR_HOOK('EPRFI2BAD_MOD:EPRFI2BAD',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE EPRFI2BAD
END MODULE EPRFI2BAD_MOD
