! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE PRFI1B_MOD
  CONTAINS
          SUBROUTINE PRFI1B(KFIRST,PSPEC,KFIELDS,KDIM,KFLDPTR)
  
  USE PARKIND1  ,ONLY : JPIM     ,JPRB
  
  USE TPM_GEN   ,ONLY : NOUT
  USE TPM_DIM   ,ONLY : R,R_NSMAX
  USE TPM_DISTR ,ONLY : D,D_NUMP,D_MYMS,D_NASM0
  USE TPM_FIELDS      ,ONLY : ZIA
  USE IEEE_ARITHMETIC
  
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
  
  INTEGER(KIND=JPIM),INTENT(IN)   :: KFIRST
  INTEGER(KIND=JPIM),INTENT(IN)   :: KFIELDS
  INTEGER(KIND=JPIM) :: KM,KMLOC
  REAL(KIND=JPRB)   ,INTENT(IN)   :: PSPEC(:,:)
  !REAL(KIND=JPRB)   ,INTENT(INOUT)  :: PIA(:,:,:)
  INTEGER(KIND=JPIM),INTENT(IN) :: KDIM
  INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL :: KFLDPTR(:)
  
  !     LOCAL INTEGER SCALARS
  INTEGER(KIND=JPIM) :: II, INM, IR, J, JFLD, ILCM, IOFF,IFLD
  
  !     ------------------------------------------------------------------
  
  !*       1.    EXTRACT FIELDS FROM SPECTRAL ARRAYS.
  !              --------------------------------------------------

#ifdef ACCGPU
  !$ACC DATA &
  !$ACC      PRESENT(D_NUMP,R_NSMAX,D_MYMS,D_NASM0) &
  !$ACC      COPYIN(PSPEC) &
  !$ACC      PRESENT(ZIA)
#endif
#ifdef OMPGPU
  !$OMP TARGET DATA MAP(PRESENT,ALLOC:D_NUMP,R_NSMAX,D_MYMS,D_NASM0,PSPEC)
#endif

#ifdef ACCGPU
  !$ACC DATA IF(PRESENT(KFLDPTR)) PRESENT(KFLDPTR)
#endif
#ifdef OMPGPU
  !$OMP TARGET DATA IF(PRESENT(KFLDPTR)) MAP(PRESENT,ALLOC:KFLDPTR)
#endif

 
  IF(PRESENT(KFLDPTR)) THEN
 
 
   !loop over wavenumber
#ifdef OMPGPU
   !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) PRIVATE(KM,ILCM,IFLD,IOFF,IR,II,INM)
#endif
#ifdef ACCGPU
   !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,ILCM,IFLD,IOFF,IR,II,INM) DEFAULT(NONE) &
   !$ACC& COPYIN(KFIELDS,KFIRST) &
   !$ACC& PRESENT(D_NUMP,D_MYMS,D_NASM0,R_NSMAX,KFLDPTR,ZIA)
#endif
   DO KMLOC=1,D_NUMP
      DO J=1,R_NSMAX+1
         DO JFLD=1,KFIELDS
            KM = D_MYMS(KMLOC)
            ILCM = R_NSMAX+1-KM
            IFLD = KFLDPTR(JFLD)
            IF (J .LE. ILCM) THEN
               IOFF = D_NASM0(KM)
               INM = IOFF+(ILCM-J)*2
               IR = KFIRST+2*(JFLD-1)
               II = IR+1
               ZIA(IR,J+2,KMLOC) = PSPEC(IFLD,INM  )
               ZIA(II,J+2,KMLOC) = PSPEC(IFLD,INM+1)
            END IF
         ENDDO
      ENDDO
 
      ! end loop over wavenumber
   END DO
#ifdef OMPGPU
   !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(KM,ILCM)
#endif
#ifdef ACCGPU
   !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM,ILCM) &
   !$ACC& COPYIN(KFIRST,KFIELDS) &
   !$ACC& PRESENT(D_NUMP,D_MYMS,R_NSMAX,ZIA)
#endif
   DO KMLOC=1,D_NUMP
      DO JFLD=KFIRST,2*KFIELDS+KFIRST-1
         KM = D_MYMS(KMLOC) 
         ILCM = R_NSMAX+1-KM
         ZIA(JFLD,1,KMLOC) = 0.0_JPRB
         ZIA(JFLD,2,KMLOC) = 0.0_JPRB
         ZIA(JFLD,ILCM+3,KMLOC) = 0.0_JPRB
      ENDDO 
      ! end loop over wavenumber
   END DO

  ELSE

   !loop over wavenumber

#ifdef OMPGPU
   !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) PRIVATE(KM,ILCM,IOFF,INM,IR,II)
#endif
#ifdef ACCGPU
   !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KMLOC,J,JFLD,KM,ILCM,IOFF,INM,IR,II) DEFAULT(NONE) &
   !$ACC& COPYIN(KFIRST,KFIELDS,KDIM) &
   !$ACC& PRESENT(ZIA,PSPEC,D_NUMP,D_MYMS,D_NASM0,R_NSMAX)
#endif
   DO KMLOC=1,D_NUMP
      DO J=1,R_NSMAX+1
         DO JFLD=1,KFIELDS
            KM = D_MYMS(KMLOC)
            ILCM = R_NSMAX+1-KM
            if (J .le. ILCM) then
               IOFF = D_NASM0(KM)
               INM = IOFF+(ILCM-J)*2
               IR = KFIRST+2*(JFLD-1)
               II = IR+1
               IF( INM .LT. KDIM ) THEN
               ZIA(IR,J+2,KMLOC) = PSPEC(JFLD,INM  )
               ZIA(II,J+2,KMLOC) = PSPEC(JFLD,INM+1)
               ENDIF
            end if
         ENDDO
      ENDDO
 
      ! end loop over wavenumber
   END DO

#ifdef OMPGPU
   !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(KM,ILCM)
#endif
#ifdef ACCGPU
   !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM,KMLOC,JFLD,ILCM) DEFAULT(NONE) &
   !$ACC& COPYIN(KFIELDS,KFIRST) &
   !$ACC& PRESENT(ZIA,D_NUMP,D_MYMS,R_NSMAX)
#endif
   DO KMLOC=1,D_NUMP
     DO JFLD=KFIRST,2*KFIELDS+KFIRST-1
         KM = D_MYMS(KMLOC) 
         ILCM = R_NSMAX+1-KM
         ZIA(JFLD,1,KMLOC) = 0.0_JPRB
         ZIA(JFLD,2,KMLOC) = 0.0_JPRB
         ZIA(JFLD,ILCM+3,KMLOC) = 0.0_JPRB
      ENDDO 
      ! end loop over wavenumber
   END DO
 
  END IF   

#ifdef ACCGPU
!!$ACC UPDATE HOST(ZIA)
!$ACC END DATA
#endif
#ifdef OMPGPU
   !$OMP END TARGET DATA
#endif
#ifdef ACCGPU
   !$ACC END DATA
#endif
#ifdef OMPGPU
   !$OMP END TARGET DATA
#endif


  !     ------------------------------------------------------------------
 
  END SUBROUTINE PRFI1B
END MODULE PRFI1B_MOD
