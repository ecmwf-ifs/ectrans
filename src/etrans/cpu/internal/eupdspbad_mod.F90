! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE EUPDSPBAD_MOD
CONTAINS
SUBROUTINE EUPDSPBAD(KM,KFIELD,POA,PSPEC,KFLDPTR)

!**** *EUPDSPBAD* - Update spectral arrays after direct Legendre transform

!     Purpose.
!     --------
!        To update spectral arrays for a fixed zonal wave-number
!         from values in POA.

!**   Interface.
!     ----------
!        CALL EUPDSPBAD(....)

!        Explicit arguments :  KM - zonal wavenumber
!        --------------------  KFIELD  - number of fields
!                              POA - work array
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
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!USE TPM_DIM
!USE TPM_FIELDS
!USE TPM_DISTR

USE TPMALD_DISTR    ,ONLY : DALD
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KM,KFIELD
REAL(KIND=JPRB)   ,INTENT(OUT)   :: POA(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPEC(:,:)
INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL :: KFLDPTR(:)

INTEGER(KIND=JPIM) :: II, INM, IR, JFLD, JN,IFLD
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

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

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EUPDSPBAD_MOD:EUPDSPBAD',0,ZHOOK_HANDLE)
POA(:,:) = 0.0_JPRB
  
IF(PRESENT(KFLDPTR)) THEN

  DO JFLD=1,KFIELD
    IR= 2*JFLD-1
    II=IR+1
    IFLD = KFLDPTR(JFLD)
!DIR$ IVDEP
!OCL NOVREC
    DO JN=1,DALD%NCPL2M(KM),2
      INM=DALD%NESM0(KM)+(JN-1)*2
      POA(JN,IR)   = PSPEC(IFLD,INM)
      POA(JN+1,IR) = PSPEC(IFLD,INM+1)
      POA(JN,II)   = PSPEC(IFLD,INM+2)
      POA(JN+1,II) = PSPEC(IFLD,INM+3)
      PSPEC(IFLD,INM  )= 0.0_JPRB
      PSPEC(IFLD,INM+1)= 0.0_JPRB
      PSPEC(IFLD,INM+2)= 0.0_JPRB
      PSPEC(IFLD,INM+3)= 0.0_JPRB
    ENDDO
  ENDDO

ELSE

  DO JN=1,DALD%NCPL2M(KM),2
    INM=DALD%NESM0(KM)+(JN-1)*2
!DIR$ IVDEP
!OCL NOVREC
    DO JFLD=1,KFIELD
      IR= 2*JFLD-1
      II=IR+1
      POA(JN,IR)   = PSPEC(JFLD,INM)
      POA(JN+1,IR) = PSPEC(JFLD,INM+1)
      POA(JN,II)   = PSPEC(JFLD,INM+2)
      POA(JN+1,II) = PSPEC(JFLD,INM+3)
      PSPEC(JFLD,INM  )= 0.0_JPRB
      PSPEC(JFLD,INM+1)= 0.0_JPRB
      PSPEC(JFLD,INM+2)= 0.0_JPRB
      PSPEC(JFLD,INM+3)= 0.0_JPRB
    ENDDO
  ENDDO

ENDIF
IF (LHOOK) CALL DR_HOOK('EUPDSPBAD_MOD:EUPDSPBAD',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE EUPDSPBAD
END MODULE EUPDSPBAD_MOD
