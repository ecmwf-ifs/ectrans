MODULE TPM_DIM

#include "tsmbkind.h"

IMPLICIT NONE

SAVE

TYPE DIM_TYPE
! SPECTRAL SPACE DIMENSIONS

  INTEGER_M :: NSMAX      ! Truncation order
  INTEGER_M :: NTMAX
  INTEGER_M :: NSPOLEG    ! Number of Legandre polynomials
  INTEGER_M :: NSPEC_G    ! Number of complex spectral coefficients (global)
  INTEGER_M :: NSPEC2_G   ! 2*NSPEC_G

! COLLOCATION GRID DIMENSIONS
  
  INTEGER_M :: NDGL       ! Number of rows of latitudes
  INTEGER_M :: NDLON      ! Maximum number of longitude points (near equator)
  INTEGER_M :: NDGNH      ! Number of rows in northern hemisphere

! Legendre transform dimensions
  INTEGER_M :: NLEI1      ! NSMAX+4+MOD(NSMAX+4+1,2) 
  INTEGER_M :: NLEI3      ! NDGNH+MOD(NDGNH+2,2)
  INTEGER_M :: NLED3
  INTEGER_M :: NLED4
END TYPE DIM_TYPE

TYPE(DIM_TYPE),ALLOCATABLE,TARGET :: DIM_RESOL(:)
TYPE(DIM_TYPE),POINTER     :: R

END MODULE TPM_DIM
