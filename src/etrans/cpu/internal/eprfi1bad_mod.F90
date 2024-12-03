! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE EPRFI1BAD_MOD
CONTAINS
SUBROUTINE EPRFI1BAD(KM,PIA,PSPEC,KFIELDS,KFLDPTR)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPMALD_DISTR    ,ONLY : DALD

!**** *EPRFI1BAD* - Prepare spectral fields for inverse Legendre transform

!     Purpose.
!     --------
!        To extract the spectral fields for a specific zonal wavenumber
!        and put them in an order suitable for the inverse Legendre           .
!        tranforms.The ordering is from NSMAX to KM for better conditioning.
!        Elements 1,2 and NLCM(KM)+1 are zeroed in preparation for computing
!        u,v and derivatives in spectral space.

!**   Interface.
!     ----------
!        *CALL* *EPRFI1BAD(...)*

!        Explicit arguments :  KM     - zonal wavenumber
!        ------------------    PIA    - spectral components for transform
!                              PSPEC  - spectral array
!                              KFIELDS  - number of fields

!        Implicit arguments :  None.
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
!        Original : 00-02-01 From PRFI1BAD in IFS CY22R1
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KM,KFIELDS
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPEC(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIA(:,:)
INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL :: KFLDPTR(:)

INTEGER(KIND=JPIM) :: II, INM, IR, J, JFLD, ILCM, IOFF, IFLD
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.    EXTRACT FIELDS FROM SPECTRAL ARRAYS.
!              --------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EPRFI1BAD_MOD:EPRFI1BAD',0,ZHOOK_HANDLE)
ILCM=DALD%NCPL2M(KM)

IOFF = DALD%NESM0(KM)

IF(PRESENT(KFLDPTR)) THEN
  DO JFLD=1,KFIELDS
    IR = 2*(JFLD-1)+1
    II = IR+1
    IFLD = KFLDPTR(JFLD)
    DO J=1,ILCM,2
      INM = IOFF+(J-1)*2

      PSPEC(IFLD,INM  ) = PSPEC(IFLD,INM  ) + PIA(J  ,IR)
      PSPEC(IFLD,INM+1) = PSPEC(IFLD,INM+1) + PIA(J+1,IR)
      PSPEC(IFLD,INM+2) = PSPEC(IFLD,INM+2) + PIA(J  ,II)
      PSPEC(IFLD,INM+3) = PSPEC(IFLD,INM+3) + PIA(J+1,II)

    ENDDO
  ENDDO
ELSE
  DO J=1,ILCM,2
    INM = IOFF+(J-1)*2
!DIR$ IVDEP
!OCL NOVREC
    DO JFLD=1,KFIELDS
      IR = 2*(JFLD-1)+1
      II = IR+1

      PSPEC(JFLD,INM  ) = PSPEC(JFLD,INM  ) + PIA(J  ,IR)
      PSPEC(JFLD,INM+1) = PSPEC(JFLD,INM+1) + PIA(J+1,IR)
      PSPEC(JFLD,INM+2) = PSPEC(JFLD,INM+2) + PIA(J  ,II)
      PSPEC(JFLD,INM+3) = PSPEC(JFLD,INM+3) + PIA(J+1,II)

    ENDDO
  ENDDO
ENDIF
IF (LHOOK) CALL DR_HOOK('EPRFI1BAD_MOD:EPRFI1BAD',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE EPRFI1BAD
END MODULE EPRFI1BAD_MOD
