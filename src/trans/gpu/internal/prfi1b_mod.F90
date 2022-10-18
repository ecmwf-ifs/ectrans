! (C) Copyright 2000- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE PRFI1B_MOD
  CONTAINS
  SUBROUTINE PRFI1B(PIA,PSPEC,KFIELDS,KDIM,KFLDPTR)
  
USE PARKIND1  ,ONLY : JPIM     ,JPRB
  
  use tpm_gen, only: nout
  USE TPM_DIM         ,ONLY : R,R_NSMAX
  USE TPM_DISTR       ,ONLY : D,D_NUMP,D_MYMS,D_NASM0
  use ieee_arithmetic
  !use openacc
  
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
  !        *CALL* *PRFI1B(...)*
  
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
  !        Original : 00-02-01 From PRFI1B in IFS CY22R1
  
  !     ------------------------------------------------------------------
  
  IMPLICIT NONE
  
  INTEGER(KIND=JPIM),INTENT(IN)   :: KFIELDS
  INTEGER(KIND=JPIM) :: KM,KMLOC
  REAL(KIND=JPRB)   ,INTENT(IN)   :: PSPEC(:,:)
  REAL(KIND=JPRB)   ,INTENT(INOUT)  :: PIA(:,:,:)
  INTEGER(KIND=JPIM),INTENT(IN) :: KDIM
  INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL :: KFLDPTR(:)
  
  !     LOCAL INTEGER SCALARS
  INTEGER(KIND=JPIM) :: II, INM, IR, JN, JFLD, ILCM, IASM0,IFLD
  
  
  !     ------------------------------------------------------------------
  
  !*       1.    EXTRACT FIELDS FROM SPECTRAL ARRAYS.
  !              --------------------------------------------------

  !$ACC DATA PRESENT(D_MYMS,D_NASM0,PIA,PSPEC)

  !$ACC DATA IF(PRESENT(KFLDPTR)) PRESENT(KFLDPTR)

  
  IF(PRESENT(KFLDPTR)) THEN
   PRINT *, "Not implemented"
   STOP 4
   
   !loop over wavenumber
   
   !$ACC PARALLEL LOOP COLLAPSE(3) DEFAULT(NONE) PRIVATE(KM,ILCM,IFLD,IASM0,IR,II,INM)
   DO KMLOC=1,D_NUMP
      DO JN=1,R_NSMAX+1
         DO JFLD=1,KFIELDS
            KM = D_MYMS(KMLOC)
            ILCM = R_NSMAX+1-KM
            IFLD = KFLDPTR(JFLD)
            IF (JN .LE. ILCM) THEN
               IASM0 = D_NASM0(KM)
               INM = IASM0+(ILCM-JN)*2
               IR = 2*(JFLD-1)+1
               II = IR+1
               PIA(IR,JN+2,KMLOC) = PSPEC(iFLD,INM  )
               PIA(II,JN+2,KMLOC) = PSPEC(iFLD,INM+1)
            END IF
         ENDDO
      ENDDO
   
      ! end loop over wavenumber
   END DO

   !$ACC PARALLEL LOOP DEFAULT(NONE) COLLAPSE(2) PRIVATE(KM,ILCM)
   DO KMLOC=1,D_NUMP
      DO JFLD=1,2*KFIELDS
         KM = D_MYMS(KMLOC) 
         ILCM = R_NSMAX+1-KM
         PIA(JFLD,1,KMLOC) = 0.0_JPRB
         PIA(JFLD,2,KMLOC) = 0.0_JPRB
         PIA(JFLD,ILCM+3,KMLOC) = 0.0_JPRB
      ENDDO 
      ! end loop over wavenumber
   END DO

ELSE

  !loop over wavenumber

  !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,IASM0,INM) DEFAULT(NONE)
  DO KMLOC=1,D_NUMP
    DO JN=0,R_NSMAX+3
      DO JFLD=1,KFIELDS
        KM = D_MYMS(KMLOC)

        IF (JN <= 1) THEN
          PIA(2*JFLD-1,JN+1,KMLOC) = 0.0_JPRB
          PIA(2*JFLD  ,JN+1,KMLOC) = 0.0_JPRB
        ELSEIF (JN <= R_NSMAX+2-KM) THEN
          IASM0 = D_NASM0(KM)
          INM = IASM0+((R_NSMAX+2-JN)-KM)*2
          PIA(2*JFLD-1,JN+1,KMLOC) = PSPEC(JFLD,INM  )
          PIA(2*JFLD  ,JN+1,KMLOC) = PSPEC(JFLD,INM+1)
        ELSEIF (JN <= R_NSMAX+3-KM) THEN
          PIA(2*JFLD-1,JN+1,KMLOC) = 0.0_JPRB
          PIA(2*JFLD  ,JN+1,KMLOC) = 0.0_JPRB
        ENDIF
      ENDDO
    ENDDO
  ENDDO

END IF

   !$ACC END DATA
   !$ACC END DATA

  !     ------------------------------------------------------------------
  
  END SUBROUTINE PRFI1B
  END MODULE PRFI1B_MOD
