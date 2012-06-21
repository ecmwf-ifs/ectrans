MODULE UPDSPAD_MOD
CONTAINS
SUBROUTINE UPDSPAD(KM,KF_UV,KF_SCALARS,POA1,POA2, &
 & PSPVOR,PSPDIV,PSPSCALAR,&
 & PSPSC3A,PSPSC3B,PSPSC2 , &
 & KFLDPTRUV,KFLDPTRSC)

!**** *UPDSPAD* - Update spectral arrays after direct Legendre transform

!     Purpose.
!     --------
!        To update the spectral arrays for a fixed zonal wave-number
!        from values in POA1 and POA2.

!**   Interface.
!     ----------
!        CALL UPDSPAD(...)

!        Explicit arguments : 
!        -------------------- 
!        KM - zonal wave-number
!        POA1 - spectral fields for zonal wavenumber KM (basic var.)
!        POA2 - spectral fields for zonal wavenumber KM (vor. div.)
!        PSPVOR - spectral vorticity
!        PSPDIV - spectral divergence
!        PSPSCALAR - spectral scalar variables

!        Implicit arguments :  
!        --------------------

!     Method.
!     -------

!     Externals.  UPDSPADB - basic transfer routine
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
!        Modified : 91-07-01 Philippe Courtier/Mats Hamrud - Rewrite
!                            for uv formulation
!        Modified : 94-08-02 R. El Khatib - interface to UPDSPADB
!        M.Hamrud : 94-11-01 New conf 'G' - vor,div->vor,div
!                            instead of u,v->vor,div
!        MPP Group: 95-10-01 Support for Distributed Memory version
!     ------------------------------------------------------------------

#include "tsmbkind.h"

USE TPM_DIM
USE TPM_TRANS
USE TPM_DISTR

USE UPDSPBAD_MOD

IMPLICIT NONE


!     DUMMY INTEGER SCALARS

INTEGER_M, INTENT(IN)  :: KM,KF_UV,KF_SCALARS

REAL_B , INTENT(OUT)  :: POA1(:,:)
REAL_B , INTENT(OUT)  :: POA2(:,:)
REAL_B  ,OPTIONAL, INTENT(INOUT) :: PSPVOR(:,:)
REAL_B  ,OPTIONAL, INTENT(INOUT) :: PSPDIV(:,:)
REAL_B  ,OPTIONAL, INTENT(INOUT) :: PSPSCALAR(:,:)
REAL_B  ,OPTIONAL, INTENT(INOUT) :: PSPSC2(:,:)
REAL_B  ,OPTIONAL, INTENT(INOUT) :: PSPSC3A(:,:,:)
REAL_B  ,OPTIONAL, INTENT(INOUT) :: PSPSC3B(:,:,:)
INTEGER_M,OPTIONAL,INTENT(IN)    :: KFLDPTRUV(:)
INTEGER_M,OPTIONAL,INTENT(IN)    :: KFLDPTRSC(:)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IVORS, IVORE, IDIVS, IDIVE, IST ,IEND, JN, ISE,IFLD,JFLD
INTEGER_M :: IDIM1,IDIM3,J3

!     ------------------------------------------------------------------

!*       1.    UPDATE FIELDS
!              -------------


!*       1.1      VORTICITY AND DIVERGENCE.

IST = 1
IF (KF_UV > 0) THEN
  IST = IST+4*KF_UV
  IVORS = 1
  IVORE = 2*KF_UV
  IDIVS = 2*KF_UV+1
  IDIVE = 4*KF_UV
  IF (KM == 0) THEN
    IF(PRESENT(KFLDPTRUV)) THEN
      DO JFLD=1,KF_UV
        IFLD = KFLDPTRUV(JFLD)
        PSPVOR(IFLD,D%NASM0(0)) = _ZERO_
        PSPDIV(IFLD,D%NASM0(0)) = _ZERO_
      ENDDO
      DO JN=0,R%NSMAX
        ISE = 1+JN*2+1
        DO JFLD=1,KF_UV
          IFLD = KFLDPTRUV(JFLD)
          PSPDIV(IFLD,ISE) = _ZERO_
          PSPVOR(IFLD,ISE) = _ZERO_
        ENDDO
      ENDDO
    ELSE
      PSPVOR(:,D%NASM0(0)) = _ZERO_
      PSPDIV(:,D%NASM0(0)) = _ZERO_
      DO JN=0,R%NSMAX
        ISE = 1+JN*2+1
        PSPDIV(:,ISE) = _ZERO_
        PSPVOR(:,ISE) = _ZERO_
      ENDDO
    ENDIF
  ENDIF
  CALL UPDSPBAD(KM,KF_UV,POA2(:,IVORS:IVORE),PSPVOR,KFLDPTRUV)
  CALL UPDSPBAD(KM,KF_UV,POA2(:,IDIVS:IDIVE),PSPDIV,KFLDPTRUV)
ENDIF

!*       1.2   SCALARS

IF (KF_SCALARS > 0) THEN
  IF(PRESENT(PSPSCALAR)) THEN
    IEND = IST+2*KF_SCALARS-1
    CALL UPDSPBAD(KM,KF_SCALARS,POA1(:,IST:IEND),PSPSCALAR,KFLDPTRSC)
  ELSE
    IF(PRESENT(PSPSC2) .AND. NF_SC2 > 0) THEN
      IDIM1 = NF_SC2
      IEND  = IST+2*IDIM1-1
      CALL UPDSPBAD(KM,IDIM1,POA1(:,IST:IEND),PSPSC2)
      IST=IST+2*IDIM1
    ENDIF
    IF(PRESENT(PSPSC3A) .AND. NF_SC3A > 0) THEN
      IDIM1=NF_SC3A
      IDIM3=UBOUND(PSPSC3A,3)
      DO J3=1,IDIM3
        IEND = IST+2*IDIM1-1
        CALL UPDSPBAD(KM,IDIM1,POA1(:,IST:IEND),PSPSC3A(:,:,J3))
        IST=IST+2*IDIM1
      ENDDO
    ENDIF
    IF(PRESENT(PSPSC3B) .AND. NF_SC3B > 0) THEN
      IDIM1=NF_SC3B
      IDIM3=UBOUND(PSPSC3B,3)
      DO J3=1,IDIM3
        IEND = IST+2*IDIM1-1
        CALL UPDSPBAD(KM,IDIM1,POA1(:,IST:IEND),PSPSC3B(:,:,J3))
        IST=IST+2*IDIM1
      ENDDO
    ENDIF
  ENDIF
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE UPDSPAD
END MODULE UPDSPAD_MOD
