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

MODULE TPM_TRANS

! Module to contain variables "local" to a specific call to a transform

!
USE PARKIND_ECTRANS,       ONLY: JPIM, JPRBT
USE GROWING_ALLOCATOR_MOD, ONLY: GROWING_ALLOCATION_TYPE

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
LOGICAL   :: LATLON     ! lat-lon output requested

!INTEGER_M :: NLEI2 ! 8*NF_UV + 2*NF_SCALARS + 2*NF_SCDERS (dimension in
                   ! inverse  Legendre transform)
!INTEGER_M :: NLED2 ! 2*NF_FS (dimension in direct Legendre transform)

!INTEGER_M :: NF_FS    ! Total number of fields in Fourier space

!INTEGER_M :: NF_GP        ! Total number of field in grid-point space
!INTEGER_M :: NF_UV_G      ! Global version of NF_UV (grid-point space)
!INTEGER_M :: NF_SCALARS_G ! Global version of NF_SCALARS (grid-point space)

INTEGER(KIND=JPIM) :: NPROMA  ! Blocking factor for gridpoint input/output
INTEGER(KIND=JPIM) :: NGPBLKS ! Number of NPROMA blocks

LOGICAL :: LGPNORM = .FALSE.  ! indicates whether transform is being done for gpnorm

! This is used in fourier space and in spectral space. It's reused among
! the transforms because we cannot reallocate - the captured graphs
! should not be modified. Hence, we keep it if it is large enough, otherwise
! we adapt the size. After 2 iterations this should lead to constant runtimes
! (the first iteration is used to get the max buffer size, the second iteration
! is going to recreate the graphs if needed)
TYPE(GROWING_ALLOCATION_TYPE) :: GROWING_ALLOCATION

END MODULE TPM_TRANS
