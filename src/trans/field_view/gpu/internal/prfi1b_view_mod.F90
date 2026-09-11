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

MODULE PRFI1B_VIEW_MOD
  CONTAINS
  SUBROUTINE PRFI1B_VIEW(PIA,YDSP)

  USE PARKIND1,        ONLY: JPIM, JPRB
  USE TPM_DIM,         ONLY: R
  USE TPM_DISTR,       ONLY: D
  USE ABORT_TRANS_MOD, ONLY: ABORT_TRANS
  USE ECTRANS_FIELD_VIEW_INTERNAL_UTIL_MOD, ONLY: SPEC_VIEW

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
  !        *CALL* *PRFI1B_VIEW(...)*

  !        Explicit arguments :  KM     - zonal wavenumber
  !        ------------------    PIA    - spectral components for transform
  !                              YDSP    - spectral arrays


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
  !        Original : 00-02-01 From PRFI1B_VIEW in IFS CY22R1

  !     ------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER(KIND=JPIM) :: KM,KMLOC
  TYPE(SPEC_VIEW), INTENT(IN) :: YDSP(:)
  REAL(KIND=JPRB)   ,INTENT(INOUT)  :: PIA(:,:,:)

  !     LOCAL INTEGER SCALARS
  INTEGER(KIND=JPIM) :: INM, IR, JN, JFLD, IASM0, IFIELDS

  !     ------------------------------------------------------------------

  !*       1.    EXTRACT FIELDS FROM SPECTRAL ARRAYS.
  !              --------------------------------------------------

  ASSOCIATE(D_NUMP=>D%NUMP, D_MYMS=>D%MYMS, D_NASM0=>D%NASM0, R_NSMAX=>R%NSMAX)
IFIELDS = SIZE(YDSP)
#ifdef ACCGPU
  !$ACC DATA PRESENT(D,D_NUMP,R,R_NSMAX,D_MYMS,D_NASM0,PIA) COPYIN(YDSP) ASYNC(1)
#endif
#ifdef OMPGPU
  !$OMP TARGET DATA MAP(PRESENT,ALLOC:D,D_NUMP,R,R_NSMAX,D_MYMS,D_NASM0,PIA,YDSP)
#endif


    !loop over wavenumber

#ifdef OMPGPU
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) DEFAULT(NONE) &
    !$OMP& PRIVATE(KM,IASM0,INM) SHARED(IFIELDS,D,R,PIA,YDSP) MAP(TO:IFIELDS)
#endif
#ifdef ACCGPU
    !$ACC PARALLEL LOOP DEFAULT(NONE) COLLAPSE(3) PRIVATE(KM,IASM0,INM) &
    !$ACC FIRSTPRIVATE(IFIELDS) &
#ifndef _CRAYFTN
    !$ACC& ASYNC(1)
#else
    !$ACC&
#endif
#endif
  DO KMLOC=1,D_NUMP
    DO JN=0,R_NSMAX+3
      DO JFLD=1,IFIELDS
        KM = D_MYMS(KMLOC)
        IF (JN <= 1) THEN
          PIA(2*JFLD-1,JN+1,KMLOC) = 0.0_JPRB
          PIA(2*JFLD  ,JN+1,KMLOC) = 0.0_JPRB
        ELSEIF (JN <= R_NSMAX+2-KM) THEN
          IASM0 = D_NASM0(KM)
          INM = IASM0+((R_NSMAX+2-JN)-KM)*2
          PIA(2*JFLD-1,JN+1,KMLOC) = YDSP(JFLD)%P(INM)
          PIA(2*JFLD  ,JN+1,KMLOC) = YDSP(JFLD)%P(INM+1)
        ELSEIF (JN <= R_NSMAX+3-KM) THEN
          PIA(2*JFLD-1,JN+1,KMLOC) = 0.0_JPRB
          PIA(2*JFLD  ,JN+1,KMLOC) = 0.0_JPRB
        ENDIF
      ENDDO
    ENDDO
  ENDDO



#ifdef ACCGPU
  !$ACC END DATA
#endif
#ifdef OMPGPU
  !$OMP END TARGET DATA
#endif

  END ASSOCIATE

  !     ------------------------------------------------------------------

  END SUBROUTINE PRFI1B_VIEW
END MODULE PRFI1B_VIEW_MOD
