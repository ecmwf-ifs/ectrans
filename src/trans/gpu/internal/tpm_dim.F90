MODULE TPM_DIM

! Module for dimensions.

USE PARKIND1  ,ONLY : JPIM     ,JPRBT

IMPLICIT NONE

SAVE

TYPE DIM_TYPE
! SPECTRAL SPACE DIMENSIONS

  INTEGER(KIND=JPIM) :: NSMAX      ! Truncation order
  INTEGER(KIND=JPIM) :: NTMAX      ! Truncation order for tendencies
  INTEGER(KIND=JPIM) :: NSPOLEG    ! Number of Legandre polynomials
  INTEGER(KIND=JPIM) :: NSPEC_G    ! Number of complex spectral coefficients (global)
  INTEGER(KIND=JPIM) :: NSPEC2_G   ! 2*NSPEC_G

! COLLOCATION GRID DIMENSIONS
  
  INTEGER(KIND=JPIM) :: NDGL       ! Number of rows of latitudes
  INTEGER(KIND=JPIM) :: NDLON      ! Maximum number of longitude points (near equator)
  INTEGER(KIND=JPIM) :: NDGNH      ! Number of rows in northern hemisphere

! Legendre transform dimensions
  INTEGER(KIND=JPIM) :: NLEI1      ! R%NSMAX+4+MOD(R%NSMAX+4+1,2)
  INTEGER(KIND=JPIM) :: NLEI3      ! R%NDGNH+MOD(R%NDGNH+2,2)
  INTEGER(KIND=JPIM) :: NLED3      ! R%NTMAX+2+MOD(R%NTMAX+3,2)
  INTEGER(KIND=JPIM) :: NLED4      ! R%NTMAX+3+MOD(R%NTMAX+4,2)

! Width of E'-zone
  INTEGER(KIND=JPIM) :: NNOEXTZL   ! Longitude direction
  INTEGER(KIND=JPIM) :: NNOEXTZG   ! Latitude direction

END TYPE DIM_TYPE

TYPE(DIM_TYPE),ALLOCATABLE,TARGET :: DIM_RESOL(:)
TYPE(DIM_TYPE),POINTER     :: R

! flat copies of above
INTEGER(KIND=JPIM) :: R_NSMAX      ! Truncation order
INTEGER(KIND=JPIM) :: R_NTMAX      ! Truncation order for tendencies
INTEGER(KIND=JPIM) :: R_NDGNH      ! Number of rows in northern hemisphere
INTEGER(KIND=JPIM) :: R_NDGL       ! Number of rows of latitudes
INTEGER(KIND=JPIM) :: R_NNOEXTZL   ! Longitude direction

END MODULE TPM_DIM
