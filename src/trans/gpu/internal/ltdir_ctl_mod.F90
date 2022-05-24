! (C) Copyright 2000- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LTDIR_CTL_MOD
  CONTAINS
  SUBROUTINE LTDIR_CTL(KF_FS,KF_UV,KF_SCALARS,FOUBUF_IN,&
   & PSPVOR,PSPDIV,PSPSCALAR, &
   & PSPSC3A,PSPSC3B,PSPSC2, &
   & KFLDPTRUV,KFLDPTRSC)
  
  !**** *LTDIR_CTL* - Control routine for direct Legendre transform
  
  !     Purpose.
  !     --------
  !        Direct Legendre transform
  
  !**   Interface.
  !     ----------
  !     CALL LTDIR_CTL(...)
  
  !     Explicit arguments :
  !     --------------------
  !     KF_FS      - number of fields in Fourier space
  !     KF_UV      - local number of spectral u-v fields
  !     KF_SCALARS - local number of scalar spectral fields
  !     PSPVOR(:,:) - spectral vorticity (output)
  !     PSPDIV(:,:) - spectral divergence (output)
  !     PSPSCALAR(:,:) - spectral scalarvalued fields (output)
  !     KFLDPTRUV(:) - field pointer for vorticity and divergence (input)
  !     KFLDPTRSC(:) - field pointer for scalarvalued fields (input)
  
  !     ------------------------------------------------------------------
  
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE PARKIND_ECTRANS  ,ONLY : JPRBT
  
  USE TPM_GEN, only: nout
  USE TPM_DIM         ,ONLY : R
  USE TPM_DISTR       ,ONLY : D
  USE TPM_GEOMETRY    ,ONLY : G
  USE TPM_FIELDS      ,ONLY : F
  
  
  USE LTDIR_MOD       ,ONLY : LTDIR
  USE TRLTOM_MOD      ,ONLY : TRLTOM, TRLTOM_CUDAAWARE
 
  USE TPM_FIELDS      ,ONLY : ZEPSNM
  USE TPM_STATS, ONLY : GSTATS => GSTATS_NVTX
  
  IMPLICIT NONE
  
  INTEGER(KIND=JPIM),INTENT(IN) :: KF_FS,KF_UV,KF_SCALARS
  REAL(KIND=JPRBT),  ALLOCATABLE, INTENT(INOUT) :: FOUBUF_IN(:)
  REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
  REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
  REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)
  REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPSC3A(:,:,:)
  REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPSC3B(:,:,:)
  REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPSC2(:,:)
  INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KFLDPTRUV(:)
  INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KFLDPTRSC(:)

  REAL(KIND=JPRBT),  ALLOCATABLE :: FOUBUF(:)
  
  INTEGER(KIND=JPIM) :: JM,IM

  CALL GSTATS(153,0)
#ifdef USE_CUDA_AWARE_MPI_FT
  WRITE(NOUT,*) 'ltdir_ctl:TRLTOM_CUDAAWARE'
  CALL TRLTOM_CUDAAWARE(FOUBUF_IN,FOUBUF,KF_FS)
#else
  CALL TRLTOM(FOUBUF_IN,FOUBUF,KF_FS)
#endif
  CALL GSTATS(153,1)
  
  ! Direct Legendre transform
  
  CALL GSTATS(103,0)
  CALL GSTATS(1645,0)

  IF (KF_FS > 0) THEN
    CALL LTDIR(FOUBUF,KF_FS,KF_UV,KF_SCALARS, &
          & PSPVOR,PSPDIV,PSPSCALAR,&
          & PSPSC3A,PSPSC3B,PSPSC2 , &
          & KFLDPTRUV,KFLDPTRSC)
  ENDIF
  CALL GSTATS(1645,1)
  
  CALL GSTATS(103,1)
  
  !     -----------------------------------------------------------------
  
  END SUBROUTINE LTDIR_CTL
  END MODULE LTDIR_CTL_MOD
