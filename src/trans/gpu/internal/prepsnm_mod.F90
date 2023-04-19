! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE PREPSNM_MOD
  CONTAINS
  SUBROUTINE PREPSNM
  
  
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
 
  USE PARKIND_ECTRANS ,ONLY : JPIM     ,JPRBT
 
  USE TPM_DIM         ,ONLY : R
  USE TPM_FIELDS      ,ONLY : F, ZEPSNM
  USE TPM_DISTR       ,ONLY : D
  USE TPM_GEN         ,ONLY : NOUT
  !
 
  IMPLICIT NONE
  
  INTEGER(KIND=JPIM)  :: KM,KMLOC
  !!REAL(KIND=JPRB),    INTENT(INOUT) :: PEPSNM(:,:)
 
  !     LOCAL INTEGER SCALARS
  INTEGER(KIND=JPIM) :: JN
  INTEGER(KIND=JPIM) :: R_NTMAX
 
 
  !     ------------------------------------------------------------------
 
  !*       1.       COPY REPSNM.
  !                 ------------
 
 
 
 
  !!!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
  !!!$ACC parallel loop
  DO KMLOC=1,D%NUMP
     KM = D%MYMS(KMLOC)
 
     IF (KM > 0) THEN
!#ifdef ACCGPU
!        !$ACC loop
!#endif
        DO JN=0,KM-1
           ZEPSNM(KMLOC,JN) = 0.0_JPRBT
        ENDDO
     ENDIF
 
     DO JN=KM,R%NTMAX+2
        ZEPSNM(KMLOC,JN) = F%REPSNM(D%NPMT(KM)+KMLOC-KM+JN)
     ENDDO
     ! end loop over wavenumber
  END DO
  !!!!$OMP END TARGET DATA
  !!!!$ACC end data
 
  !     ------------------------------------------------------------------
 
  END SUBROUTINE PREPSNM
 
 END MODULE PREPSNM_MOD  
