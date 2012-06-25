MODULE TPM_FFT
USE PARKIND1  ,ONLY : JPIM     ,JPRB

! Module for Fourier transforms.

IMPLICIT NONE

SAVE

TYPE FFT_TYPE
  REAL(KIND=JPRB)   ,POINTER :: TRIGS(:,:) ! list of trigonometric function values
  INTEGER(KIND=JPIM),POINTER :: NFAX(:,:)  ! list of factors of truncation
END TYPE FFT_TYPE

TYPE(FFT_TYPE),ALLOCATABLE,TARGET :: FFT_RESOL(:)
TYPE(FFT_TYPE),POINTER     :: T

END MODULE TPM_FFT
