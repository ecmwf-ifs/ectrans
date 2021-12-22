MODULE PREPSNM_MOD
  CONTAINS
  SUBROUTINE PREPSNM(PEPSNM)
  
  
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
  
  USE PARKIND1  ,ONLY : JPIM     ,JPRBT
  
  USE TPM_DIM         ,ONLY : R
  USE TPM_FIELDS      ,ONLY : F
  USE TPM_DISTR       ,ONLY : D
  use tpm_gen, only: nout
  !
  
  IMPLICIT NONE
  
  INTEGER(KIND=JPIM)  :: KM,KMLOC
  REAL(KIND=JPRBT),    INTENT(OUT) :: PEPSNM(1:d%nump,0:R%NTMAX+2)
  !REAL(KIND=JPRBT),    INTENT(OUT) :: PEPSNM(:,:)
  
  !     LOCAL INTEGER SCALARS
  INTEGER(KIND=JPIM) :: JN
  INTEGER(KIND=JPIM) :: R_NTMAX
  
  
  !     ------------------------------------------------------------------
  
  !*       1.       COPY REPSNM.
  !                 ------------
  
  
  
  R_NTMAX = R%NTMAX
  
  !$ACC data &
  !$ACC& COPYIN(D,F,D%NPMT,D%NUMP,D%MYMS,F%REPSNM) &
  !$ACC& present(PEPSNM)
  
  !$ACC parallel loop
  DO KMLOC=1,D%NUMP
     KM = D%MYMS(KMLOC)
     
     
     IF (KM > 0) THEN
        !PEPSNM(0:KM-1) = 0.0_JPRBT
        !$ACC loop
        DO JN=0,KM-1
           PEPSNM(KMLOC,JN) = 0.0_JPRBT
        ENDDO
     ENDIF
  
     !$ACC loop
     DO JN=KM,R_NTMAX+2
        PEPSNM(KMLOC,JN) = F%REPSNM(D%NPMT(KM)+KMLOC-KM+JN)
     ENDDO
     ! end loop over wavenumber
  END DO

  !$ACC end data
  
  !     ------------------------------------------------------------------
  
  END SUBROUTINE PREPSNM
  
  END MODULE PREPSNM_MOD  
