! (C) Copyright 2000- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LTINV_CTL_MOD
  CONTAINS
  SUBROUTINE LTINV_CTL(KF_UV,KF_SCALARS,&
   & PSPVOR,PSPDIV,PSPSCALAR,PSPSC3A,PSPSC3B,PSPSC2,FOUBUF,KFIELD,&
   & KFLDPTRUV,KFLDPTRSC)
  
  !**** *LTINV_CTL* - Control routine for inverse Legandre transform.
  
  !     Purpose.
  !     --------
  !        Control routine for the inverse LEGENDRE transform
  
  !**   Interface.
  !     ----------
  !     CALL INV_TRANS_CTL(...)
  !     KF_UV        - local number of spectral u-v fields
  !     KF_SCALARS   - local number of scalar spectral fields
  !     PSPVOR(:,:)  - spectral vorticity (input)
  !     PSPDIV(:,:)  - spectral divergence (input)
  !     PSPSCALAR(:,:) - spectral scalarvalued fields (input)
  !     KFLDPTRUV(:) - field pointer array for vor./div.
  !     KFLDPTRSC(:) - field pointer array for PSPSCALAR
  
  !     Method.
  !     -------
  
  !     Externals.
  !     ----------
  !
  
  !     Author.
  !     -------
  !        Mats Hamrud *ECMWF*
  
  !     Modifications.
  !     --------------
  !        Original : 00-06-03
  
  !     ------------------------------------------------------------------
  
  USE PARKIND1  ,ONLY : JPIM     ,JPRB
  
  USE TPM_GEN, only: nout
  USE TPM_DIM         ,ONLY : R
  USE TPM_DISTR       ,ONLY : D
  USE TPM_GEOMETRY    ,ONLY : G
  
  USE TPM_FLT
  
  USE LTINV_MOD       ,ONLY : LTINV
  USE TRMTOL_MOD      ,ONLY : TRMTOL, TRMTOL_CUDAAWARE
  USE TPM_STATS, ONLY : GSTATS => GSTATS_NVTX

  IMPLICIT NONE

  INTEGER(KIND=JPIM),INTENT(IN) :: KF_UV,KF_SCALARS
  REAL(KIND=JPRB) ,OPTIONAL, INTENT(IN)  :: PSPVOR(:,:)
  REAL(KIND=JPRB) ,OPTIONAL, INTENT(IN)  :: PSPDIV(:,:)
  REAL(KIND=JPRB) ,OPTIONAL, INTENT(IN)  :: PSPSCALAR(:,:)
  REAL(KIND=JPRB) ,OPTIONAL, INTENT(IN)  :: PSPSC3A(:,:,:)
  REAL(KIND=JPRB) ,OPTIONAL, INTENT(IN)  :: PSPSC3B(:,:,:)
  REAL(KIND=JPRB) ,OPTIONAL, INTENT(IN)  :: PSPSC2(:,:)
  REAL(KIND=JPRB), ALLOCATABLE, INTENT(OUT)  :: FOUBUF(:)
  INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KFLDPTRUV(:)
  INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KFLDPTRSC(:)
  INTEGER(KIND=JPIM),INTENT(OUT) :: KFIELD
  REAL(KIND=JPRB), ALLOCATABLE :: FOUBUF_IN(:)

  
  !     ------------------------------------------------------------------
  
  END SUBROUTINE LTINV_CTL
  END MODULE LTINV_CTL_MOD
