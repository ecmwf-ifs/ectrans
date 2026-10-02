! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE PRFI1B_VIEW_MOD
CONTAINS
SUBROUTINE PRFI1B_VIEW(KM,PIA,YDSP)

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DIM         ,ONLY : R
USE TPM_DISTR       ,ONLY : D
USE ECTRANS_FIELD_VIEW_INTERNAL_UTIL_MOD, ONLY: SPEC_VIEW

!**** *PRFI1* - Prepare spectral fields for inverse Legendre transform

!     Purpose.
!     --------
!        To extract the spectral fields for a specific zonal wavenumber
!        and put them in an order suitable for the inverse Legendre           .
!        tranforms.The ordering is from NSMAX to KM for better conditioning.
!        Elements 1,2 and NLCM(KM)+1 are zeroed in preparation for computing
!        u,v and derivatives in spectral space.

!**   Interface.
!     ----------
!        *CALL* *PRFI1B_VIEW(...)*

!        Explicit arguments :  KM     - zonal wavenumber
!        ------------------    PIA    - spectral components for transform
!                              YDSP    - spectral arrays



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
!        Original : 00-02-01 From PRFI1B_VIEW in IFS CY22R1

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)   :: KM
TYPE(SPEC_VIEW), INTENT(IN) :: YDSP(:)
REAL(KIND=JPRB)   ,INTENT(OUT)  :: PIA(:,:)


!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: II, INM, IR, J, JFLD, ILCM, IOFF,IFLD,IFIELDS


!     ------------------------------------------------------------------

!*       1.    EXTRACT FIELDS FROM SPECTRAL ARRAYS.
!              --------------------------------------------------

ILCM = R%NSMAX+1-KM
IOFF = D%NASM0(KM)
IFIELDS = SIZE(YDSP)

DO J=1,ILCM
  INM = IOFF+(ILCM-J)*2
  !DIR$ IVDEP
  !OCL NOVREC
  DO JFLD=1,IFIELDS
    IR = 2*(JFLD-1)+1
    II = IR+1
    PIA(J+2,IR) = YDSP(JFLD)%P(INM)
    PIA(J+2,II) = YDSP(JFLD)%P(INM+1)
  ENDDO
ENDDO

DO JFLD=1,2*IFIELDS
  PIA(1,JFLD) = 0.0_JPRB
  PIA(2,JFLD) = 0.0_JPRB
  PIA(ILCM+3:,JFLD) = 0.0_JPRB
ENDDO


!     ------------------------------------------------------------------

END SUBROUTINE PRFI1B_VIEW
END MODULE PRFI1B_VIEW_MOD
