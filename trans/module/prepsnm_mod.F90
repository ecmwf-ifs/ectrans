MODULE PREPSNM_MOD
CONTAINS
SUBROUTINE PREPSNM(KM,KMLOC,PEPSNM)

#ifdef DOC

!**** *PREPSNM* - Prepare REPSNM for wavenumber KM

!     Purpose.
!     --------
!        Copy the REPSNM values for specific zonal wavenumber M
!        to work array

!**   Interface.
!     ----------
!        CALL PREPSNM(...)

!        Explicit arguments :  KM - zonal wavenumber
!        -------------------   KMLOC - local zonal wavenumber
!                              PEPSNM - REPSNM for zonal
!                                      wavenumber KM

!        Implicit arguments :  
!        --------------------

!     Method.
!     -------


!     Reference.
!     ----------


!     Author.
!     -------
!        Lars Isaksen *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From LTINV in IFS CY22R1

!     ------------------------------------------------------------------
#endif

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DIM
USE TPM_FIELDS
USE TPM_DISTR

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)  :: KM,KMLOC
REAL(KIND=JPRB),    INTENT(OUT) :: PEPSNM(0:R%NTMAX+2)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: JN

!     ------------------------------------------------------------------

!*       1.       COPY REPSNM.
!                 ------------


IF (KM > 0) THEN
  PEPSNM(0:KM-1) = 0.0_JPRB
ENDIF

DO JN=KM,R%NTMAX+2
  PEPSNM(JN) = F%REPSNM(D%NPMT(KM)+KMLOC-KM+JN)
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE PREPSNM
END MODULE PREPSNM_MOD

