MODULE LEINV_MOD
CONTAINS
SUBROUTINE LEINV(KM,KFC,KF_OUT_LT,PIA,PAOA1,PSOA1)

!**** *LEINV* - Inverse Legendre transform.

!     Purpose.
!     --------
!        Inverse Legendre tranform of all variables(kernel).

!**   Interface.
!     ----------
!        CALL LEINV(...)

!        Explicit arguments :  KM - zonal wavenumber (input-c)
!        --------------------  KFC - number of fields to tranform (input-c)
!                              PIA - spectral fields
!                              for zonal wavenumber KM (input)
!                              PAOA1 - antisymmetric part of Fourier
!                              fields for zonal wavenumber KM (output)
!                              PSOA1 - symmetric part of Fourier
!                              fields for zonal wavenumber KM (output)

!        Implicit arguments :  None.
!        --------------------

!     Method.
!     -------

!     Externals.   MXMAOP - calls SGEMVX (matrix multiply)
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From LEINV in IFS CY22R1

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DIM
USE TPM_GEOMETRY
USE TPM_TRANS
USE TPM_FIELDS
USE TPM_DISTR

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)  :: KM
INTEGER(KIND=JPIM), INTENT(IN)  :: KFC
INTEGER(KIND=JPIM), INTENT(IN)  :: KF_OUT_LT
REAL(KIND=JPRB),    INTENT(IN)  :: PIA(:,:)
REAL(KIND=JPRB),    INTENT(OUT) :: PSOA1(:,:)
REAL(KIND=JPRB),    INTENT(OUT) :: PAOA1(:,:)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: IA, IDGLU, IFC, ILA, ILS, IS, ISKIP, ISL, J1, JGL,IOAD1


!     ------------------------------------------------------------------

!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------

!*       1.1      PREPARATIONS.

ISL = MAX(R%NDGNH-G%NDGLU(KM)+1,1)

IA  = 1+MOD(R%NSMAX-KM+2,2)
IS  = 1+MOD(R%NSMAX-KM+1,2)
ILA = (R%NSMAX-KM+2)/2
ILS = (R%NSMAX-KM+3)/2
ISL = MAX(R%NDGNH-G%NDGLU(KM)+1,1)
IOAD1 = 2*KF_OUT_LT

IF(KM == 0)THEN
  ISKIP = 2
  IFC = KFC/2
  DO J1=2,KFC,2
    DO JGL=ISL,R%NDGNH
      PSOA1(J1,JGL) = 0.0_JPRB
      PAOA1(J1,JGL) = 0.0_JPRB
    ENDDO
  ENDDO
ELSE
  ISKIP = 1
  IFC = KFC
ENDIF

!*       1.2      ANTISYMMETRIC PART.

IDGLU = MIN(R%NDGNH,G%NDGLU(KM))

IF( IDGLU > 0 ) THEN

  CALL MXMAOP(F%RPNM(ISL,D%NPMS(KM)+IA),1,2*R%NLEI3,PIA(1+IA,1),2,&
              &R%NLEI1*ISKIP,PAOA1(1,ISL),IOAD1,ISKIP,&
              &IDGLU,ILA,IFC)

!*       1.3      SYMMETRIC PART.

  CALL MXMAOP(F%RPNM(ISL,D%NPMS(KM)+IS),1,2*R%NLEI3,PIA(1+IS,1),2,&
              &R%NLEI1*ISKIP,PSOA1(1,ISL),IOAD1,ISKIP,&
              &IDGLU,ILS,IFC)

ENDIF

!     ------------------------------------------------------------------


END SUBROUTINE LEINV
END MODULE LEINV_MOD
