MODULE UPDSPB_MOD
CONTAINS
SUBROUTINE UPDSPB(KM,KFIELD,POA,PSPEC)


!**** *UPDSPB* - Update spectral arrays after direct Legendre transform

!     Purpose.
!     --------
!        To update spectral arrays for a fixed zonal wave-number
!         from values in POA.

!**   Interface.
!     ----------
!        CALL UPDSPB(....)

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
!     ------------------------------------------------------------------

#include "tsmbkind.h"

USE TPM_DIM
USE TPM_FIELDS
USE TPM_DISTR

IMPLICIT NONE

INTEGER_M,INTENT(IN)  :: KM,KFIELD
REAL_B   ,INTENT(IN)  :: POA(:,:)
REAL_B   ,INTENT(OUT) :: PSPEC(:,:)

!     LOCAL INTEGER SCALARS
INTEGER_M :: II, INM, IR, JFLD, JN, ISMAX, ITMAX, IASM0


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

!*       1.1   KM=0

IF(KM == 0) THEN
  DO JN=ITMAX+2-ISMAX,ITMAX+2-KM
    INM = IASM0+(ITMAX+2-JN)*2
!DIR$ IVDEP
!OCL NOVREC
    DO JFLD=1,KFIELD
      IR = 2*JFLD-1
      PSPEC(JFLD,INM)   = POA(JN,IR)
      PSPEC(JFLD,INM+1) = _ZERO_
    ENDDO
  ENDDO

!*       1.2   KM!=0

ELSE
  DO JN=ITMAX+2-ISMAX,ITMAX+2-KM
    INM = IASM0+((ITMAX+2-JN)-KM)*2
!DIR$ IVDEP
!OCL NOVREC
    DO JFLD=1,KFIELD
      IR = 2*JFLD-1
      II = IR+1
      PSPEC(JFLD,INM)   = POA(JN,IR)
      PSPEC(JFLD,INM+1) = POA(JN,II)
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE UPDSPB
END MODULE UPDSPB_MOD
