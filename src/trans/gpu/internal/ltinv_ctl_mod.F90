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
  SUBROUTINE LTINV_CTL(KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS,&
   & PSPVOR,PSPDIV,PSPSCALAR,&
   & PSPSC3A,PSPSC3B,PSPSC2,&
   & KFLDPTRUV,KFLDPTRSC,FSPGL_PROC)
  
  !**** *LTINV_CTL* - Control routine for inverse Legandre transform.
  
  !     Purpose.
  !     --------
  !        Control routine for the inverse LEGENDRE transform
  
  !**   Interface.
  !     ----------
  !     CALL INV_TRANS_CTL(...)
  !     KF_OUT_LT    - number of fields coming out from inverse LT
  !     KF_UV        - local number of spectral u-v fields
  !     KF_SCALARS   - local number of scalar spectral fields
  !     KF_SCDERS    - local number of derivatives of scalar spectral fields
  !     PSPVOR(:,:)  - spectral vorticity (input)
  !     PSPDIV(:,:)  - spectral divergence (input)
  !     PSPSCALAR(:,:) - spectral scalarvalued fields (input)
  !     KFLDPTRUV(:) - field pointer array for vor./div.
  !     KFLDPTRSC(:) - field pointer array for PSPSCALAR
  !     FSPGL_PROC  - external procedure to be executed in fourier space
  !                   before transposition
  
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
  USE TPM_TRANS       ,ONLY : FOUBUF, FOUBUF_IN
  USE TPM_DISTR       ,ONLY : D
  USE TPM_GEOMETRY    ,ONLY : G
  USE TPM_FIELDS      ,ONLY : ZGTF_START, ZGTF_START_INDEX_NSDERS
  
  USE TPM_FLT
  
  USE LTINV_MOD       ,ONLY : LTINV
  USE TRMTOL_MOD      ,ONLY : TRMTOL, TRMTOL_CUDAAWARE
  
  IMPLICIT NONE
  
  INTEGER(KIND=JPIM),INTENT(IN) :: KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS
  REAL(KIND=JPRB) ,OPTIONAL, INTENT(IN)  :: PSPVOR(:,:)
  REAL(KIND=JPRB) ,OPTIONAL, INTENT(IN)  :: PSPDIV(:,:)
  REAL(KIND=JPRB) ,OPTIONAL, INTENT(IN)  :: PSPSCALAR(:,:)
  REAL(KIND=JPRB) ,OPTIONAL, INTENT(IN)  :: PSPSC3A(:,:,:)
  REAL(KIND=JPRB) ,OPTIONAL, INTENT(IN)  :: PSPSC3B(:,:,:)
  REAL(KIND=JPRB) ,OPTIONAL, INTENT(IN)  :: PSPSC2(:,:)
  INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KFLDPTRUV(:)
  INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KFLDPTRSC(:)
  EXTERNAL  FSPGL_PROC
  OPTIONAL  FSPGL_PROC
  
  INTEGER(KIND=JPIM) :: JM,IM,IBLEN,ILEI2,IDIM1, i, j
  
  
  !$ACC DATA CREATE(FOUBUF_IN) PRESENT(FOUBUF)
  
  CALL GSTATS(102,0)
  ILEI2 = 8*KF_UV + 2*KF_SCALARS + 2*KF_SCDERS
  IDIM1 = 2*KF_OUT_LT
  IBLEN = D%NLENGT0B*2*KF_OUT_LT

  IF(KF_OUT_LT > 0) THEN
    CALL GSTATS(1647,0)

      ! from PSPXXX to FOUBUF_IN
    CALL LTINV(KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS,ILEI2,IDIM1,&
        & PSPVOR,PSPDIV,PSPSCALAR ,&
        & PSPSC3A,PSPSC3B,PSPSC2 , &
        & KFLDPTRUV,KFLDPTRSC,FSPGL_PROC)

    CALL GSTATS(1647,1)
  ENDIF
  CALL GSTATS(102,1)

  CALL GSTATS(152,0)
  ! from FOUBUF_IN to FOUBUF
#ifdef USE_CUDA_AWARE_MPI_FT
  WRITE(NOUT,*) 'ltinv_ctl:TRMTOL_CUDAAWARE'
  CALL TRMTOL_CUDAAWARE(FOUBUF_IN,FOUBUF,ZGTF_START(ZGTF_START_INDEX_NSDERS+1)-1)
#else
  CALL TRMTOL(FOUBUF_IN,FOUBUF,ZGTF_START(ZGTF_START_INDEX_NSDERS+1)-1)
#endif
  CALL GSTATS(152,1)
  !$ACC END DATA
  
  !     ------------------------------------------------------------------
  
  END SUBROUTINE LTINV_CTL
  END MODULE LTINV_CTL_MOD
