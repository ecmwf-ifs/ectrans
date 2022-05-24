! (C) Copyright 2000- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE SUFFT_MOD
  CONTAINS
  SUBROUTINE SUFFT
  
  USE PARKIND1  ,ONLY : JPIM
  
  USE TPM_DIM         ,ONLY : R
  USE TPM_GEN         ,ONLY : NOUT, NPRINTLEV
  USE TPM_DISTR       ,ONLY : D, MYSETW
  USE TPM_GEOMETRY    ,ONLY : G
  USE TPM_FFT         ,ONLY : T
#ifdef WITH_FFTW
  USE TPM_FFTW        ,ONLY : TW, INIT_PLANS_FFTW
#endif
  !
  
  IMPLICIT NONE
  
  INTEGER(KIND=JPIM) :: JGL,IGLG
  LOGICAL :: LLP1,LLP2
  
  !     ------------------------------------------------------------------
  
  IF(.NOT.D%LGRIDONLY) THEN
  
    LLP1 = NPRINTLEV>0
    LLP2 = NPRINTLEV>1
    IF(LLP1) WRITE(NOUT,*) '=== ENTER ROUTINE SUFFT ==='

#ifdef WITH_FFTW
    IF(TW%LFFTW)THEN
      CALL INIT_PLANS_FFTW(R%NDLON)
    ELSE
      NULLIFY(TW%FFTW_PLANS)
      !NULLIFY(TW%N_PLANS)
    ENDIF
#endif
  
  ENDIF
  
  !     ------------------------------------------------------------------
  
  9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)
  
  END SUBROUTINE SUFFT
  END MODULE SUFFT_MOD  
