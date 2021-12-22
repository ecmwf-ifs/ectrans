MODULE TPM_FFT
USE PARKIND1  ,ONLY : JPIM     ,JPRBT
USE BLUESTEIN_MOD ,ONLY : FFTB_TYPE

! Module for Fourier transforms.

IMPLICIT NONE

SAVE

TYPE FFT_TYPE
  REAL(KIND=JPRBT)   ,ALLOCATABLE :: TRIGS(:,:) ! list of trigonometric function values
  INTEGER(KIND=JPIM),ALLOCATABLE :: NFAX(:,:)  ! list of factors of truncation
  LOGICAL                        :: LBLUESTEIN=.FALSE. ! logical indicating whether some
                                                       ! latitudes require bluestein algorithm
                                                       ! with prime factors that are not 2,3,or 5
  LOGICAL,ALLOCATABLE            :: LUSEFFT992(:) ! describes which FFT algorithm to be used
                                                  ! T=use FFT992 F=use bluestein
END TYPE FFT_TYPE

TYPE(FFT_TYPE),ALLOCATABLE,TARGET :: FFT_RESOL(:)
TYPE(FFT_TYPE),POINTER     :: T

TYPE(FFTB_TYPE),ALLOCATABLE,TARGET :: FFTB_RESOL(:)
TYPE(FFTB_TYPE),POINTER    :: TB

END MODULE TPM_FFT
