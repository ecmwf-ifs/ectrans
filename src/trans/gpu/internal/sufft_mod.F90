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

MODULE SUFFT_MOD
  CONTAINS
  SUBROUTINE SUFFT

  USE PARKIND1,     ONLY: JPIM
  USE TPM_DIM,      ONLY: R
  USE TPM_GEN,      ONLY: NOUT, NPRINTLEV
  USE TPM_DISTR,    ONLY: D
  USE TPM_HICFFT,   ONLY: INIT_PLANS_FFT
  !

  IMPLICIT NONE

  INTEGER(KIND=JPIM) :: JGL,IGLG
  LOGICAL :: LLP1,LLP2

  !     ------------------------------------------------------------------

  IF(.NOT.D%LGRIDONLY) THEN

    LLP1 = NPRINTLEV>0
    LLP2 = NPRINTLEV>1
    IF(LLP1) WRITE(NOUT,*) '=== ENTER ROUTINE SUFFT ==='

    CALL INIT_PLANS_FFT(R%NDLON)

  ENDIF

  !     ------------------------------------------------------------------

  9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

  END SUBROUTINE SUFFT
END MODULE SUFFT_MOD
