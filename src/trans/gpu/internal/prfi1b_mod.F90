! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
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
  
  USE PARKIND1,        ONLY: JPIM, JPRB
  USE TPM_DIM,         ONLY: R
  USE TPM_DISTR,       ONLY: D
  USE ABORT_TRANS_MOD, ONLY: ABORT_TRANS
  
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
  INTEGER(KIND=JPIM) :: II, INM, IR, JN, JFLD, IASM0
  
  !     ------------------------------------------------------------------
  
  !*       1.    EXTRACT FIELDS FROM SPECTRAL ARRAYS.
  !              --------------------------------------------------

  ASSOCIATE(D_NUMP=>D%NUMP, D_MYMS=>D%MYMS, D_NASM0=>D%NASM0, R_NSMAX=>R%NSMAX)

#ifdef ACCGPU
  !$ACC DATA PRESENT(D,D_NUMP,R,R_NSMAX,D_MYMS,D_NASM0,PIA,PSPEC) ASYNC(1)
#endif
 
  IF(PRESENT(KFLDPTR)) THEN
 
    CALL ABORT_TRANS("KFLDPTR not implemented for GPU")

  ELSE

    !loop over wavenumber

#ifdef OMPGPU
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) DEFAULT(NONE) &
    !$OMP& PRIVATE(KM,IASM0,INM) SHARED(KFIELDS,KDIM,D,R,PIA,PSPEC) MAP(TO:KFIELDS)
#endif
#ifdef ACCGPU
    !$ACC PARALLEL LOOP DEFAULT(NONE) COLLAPSE(3) PRIVATE(KM,IASM0,INM) &
    !$ACC FIRSTPRIVATE(KFIELDS,KDIM) &
#ifndef _CRAYFTN
    !$ACC& ASYNC(1)
#else
    !$ACC&
#endif
#endif
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

  ENDIF

#ifdef ACCGPU
  !$ACC END DATA
#endif

  END ASSOCIATE

  !     ------------------------------------------------------------------

  END SUBROUTINE PRFI1B
END MODULE PRFI1B_MOD
