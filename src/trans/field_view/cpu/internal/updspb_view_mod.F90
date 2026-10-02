! (C) Copyright 1988- ECMWF.
! (C) Copyright 1988- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE UPDSPB_VIEW_MOD
CONTAINS
SUBROUTINE UPDSPB_VIEW(KM,POA,YDSPEC)


!**** *UPDSPB_VIEW* - Update spectral arrays after direct Legendre transform

!     Purpose.
!     --------
!        To update spectral arrays for a fixed zonal wave-number
!         from values in POA.

!**   Interface.
!     ----------
!        CALL UPDSPB_VIEW(....)

!        Explicit arguments :  KM - zonal wavenumber
!        --------------------  POA - work array!
!                              PSPEC - spectral array

!        Implicit arguments :  None
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 88-02-02
!        D. Giard : 93-03-19 truncations NSMAX and NTMAX (see NOTE)
!        R. El Khatib : 94-08-02 Replace number of fields by indexes of the
!                       first and last field
!        L. Isaksen : 95-06-06 Reordering of spectral arrays
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DIM         ,ONLY : R
USE TPM_DISTR       ,ONLY : D
USE ECTRANS_FIELD_VIEW_INTERNAL_UTIL_MOD, ONLY : SPEC_VIEW

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KM
REAL(KIND=JPRB)   ,INTENT(IN)  :: POA(:,:)
TYPE(SPEC_VIEW) :: YDSPEC(:)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: IR, INM, JFLD, JN, ISMAX, ITMAX, IASM0,IFLD
REAL(KIND=JPRB) :: INM1

!     ------------------------------------------------------------------

!*       0.    NOTE.
!              -----

! The following transfer reads :
! SPEC(k,NASM0(m)+NLTN(n)*2)  =POA(nn,2*k-1) (real part)
! SPEC(k,NASM0(m)+NLTN(n)*2+1)=POA(nn,2*k  ) (imaginary part)
! with n from m to NSMAX
! and nn=NTMAX+2-n from NTMAX+2-m to NTMAX+2-NSMAX.
! NLTN(m)=NTMAX+2-m : n=NLTN(nn),nn=NLTN(n)
! nn is the loop index.



!*       1.    UPDATE SPECTRAL FIELDS.
!              -----------------------
ISMAX = R%NSMAX
ITMAX = R%NTMAX
IASM0 = D%NASM0(KM)

INM1 = 0.0_JPRB
DO JN=ITMAX+2-ISMAX,ITMAX+2-KM
      INM = IASM0+((ITMAX+2-JN)-KM)*2
!DIR$ IVDEP
!OCL NOVREC
  DO JFLD=1,SIZE(YDSPEC)
    IR = 2*JFLD-1
    YDSPEC(JFLD)%P(INM) = POA(JN,IR)
    IF (KM /= 0) INM1 = POA(JN,IR+1)
    YDSPEC(JFLD)%P(INM+1) = INM1
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE UPDSPB_VIEW
END MODULE UPDSPB_VIEW_MOD
