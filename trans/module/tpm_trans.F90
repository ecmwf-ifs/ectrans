MODULE TPM_TRANS

! Module to contain variables "local" to a specific call to a transform

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!INTEGER_M :: NF_UV      ! Number of u-v fields (spectral/fourier space)
!INTEGER_M :: NF_SCALARS ! Number of scalar fields (spectral/fourier space)
!INTEGER_M :: NF_SCDERS  ! Number of fields for derivatives of scalars
                        ! (inverse transform, spectral/fourier space)
!INTEGER_M :: NF_OUT_LT  ! Number of fields that comes out of Inverse 
                        ! Legendre transform
INTEGER(KIND=JPIM) :: NF_SC2  ! Number of fields in "SPSC2" arrays.
INTEGER(KIND=JPIM) :: NF_SC3A ! Number of fields in "SPSC3A" arrays.
INTEGER(KIND=JPIM) :: NF_SC3B ! Number of fields in "SPSC3B" arrays.

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

REAL(KIND=JPRB), ALLOCATABLE :: FOUBUF_IN(:)  ! Fourier buffer 
REAL(KIND=JPRB), ALLOCATABLE :: FOUBUF(:)     ! Fourier buffer

INTEGER(KIND=JPIM) :: NPROMA  ! Blocking factor for gridpoint input/output
INTEGER(KIND=JPIM) :: NGPBLKS ! Number of NPROMA blocks

END MODULE TPM_TRANS
