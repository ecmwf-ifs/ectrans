MODULE TPM_TRANS

! Module to contain variables "local" to a specific call to a transform

#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!INTEGER_M :: NF_UV      ! Number of u-v fields (spectral/fourier space)
!INTEGER_M :: NF_SCALARS ! Number of scalar fields (spectral/fourier space)
!INTEGER_M :: NF_SCDERS  ! Number of fields for derivatives of scalars
                        ! (inverse transform, spectral/fourier space)
!INTEGER_M :: NF_OUT_LT  ! Number of fields that comes out of Inverse 
                        ! Legendre transform
INTEGER_M :: NF_SC2
INTEGER_M :: NF_SC3A
INTEGER_M :: NF_SC3B

!LOGICAL   :: LUV        ! uv fields requested
!LOGICAL   :: LSCALAR    ! scalar fields requested
LOGICAL   :: LVORGP     ! vorticity requested
LOGICAL   :: LDIVGP     ! divergence requested
LOGICAL   :: LUVDER     ! E-W derivatives of U and V requested
LOGICAL   :: LSCDERS    ! derivatives of scalar variables are req.

!INTEGER_M :: NLEI2 ! 8*NF_UV + 2*NF_SCALARS + 2*NF_SCDERS (dimension in 
                   ! inverse  Legendre transform)
!INTEGER_M :: NLED2 ! 2*NF_FS (dimension in direct Legendre transform)

!INTEGER_M :: NF_FS    ! Total number of fields in Fourier space

!INTEGER_M :: NF_GP        ! Total number of field in grid-point space
!INTEGER_M :: NF_UV_G      ! Global version of NF_UV (grid-point space)
!INTEGER_M :: NF_SCALARS_G ! Global version of NF_SCALARS (grid-point space)

REAL_B, ALLOCATABLE :: FOUBUF_IN(:)  ! Fourier buffer 
REAL_B, ALLOCATABLE :: FOUBUF(:)     ! Fourier buffer

INTEGER_M :: NPROMA  ! Blocking factor for gridpoint input/output
INTEGER_M :: NGPBLKS ! Number of NPROMA blocks

END MODULE TPM_TRANS
