! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

INTERFACE
SUBROUTINE SPECNORM(PNORM,PSPEC,KVSET,KMASTER,KRESOL,PMET)

! begin_doc_block
! ## `SPECNORM`
!
! ### Signature
!
! ```f90
! SUBROUTINE SPECNORM(PNORM, PSPEC, KVSET, KMASTER, KRESOL, PMET)
! ```
!
! ### Purpose
!
! This subroutine computes the spectral norms of the fields in the given input array and returns the
! values to a specific MPI task.
!
! ### `OPTIONAL, INTENT(IN)` arguments
!
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PSPEC(:,:)`  
!   The spectral array to compute the norm of. This argument is optional because not all MPI tasks  
!   necessarily possess a field, depending on the contents of `KVSET`.  
!   Dimensions: (number of fields on this MPI task, spectral dimension)  
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KVSET(:)`  
!   Array specifying which MPI task each of the fields belongs to.  
!   Dimensions: (number of global fields)
!   Dimensions: (number of fields across all MPI tasks)  
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KMASTER`  
!   The MPI task for whom `PNORM` will contain values.  
!   *Default*: `1`
! - `REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: PMET(:)`  
!   A multiplier applied in the total wavenumber direction as the norm is calculated.  
!   Dimensions: (`0:R%NSMAX`) (`R%NSMAX` is the same as `KSMAX` passed to `SETUP_TRANS`)  
!   *Default*: `1.0`
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KRESOL`  
!   Resolution handle returned by original call to `SETUP_TRANS`.  
!   *Default*: `1` (i.e. first resolution handle)
!
! ### `INTENT(OUT)` arguments
!
! - `REAL(KIND=JPRB), INTENT(OUT) :: PNORM(:)`  
!   An array containing the spectral norms of all input spectral fields.  
!   Note that all MPI tasks must pass this array - it is not optional. Only task `KMASTER` needs to
!   actually allocate its array though.  
!   Dimensions: (number of fields across all MPI tasks)
! end_doc_block

USE PARKIND1  ,ONLY : JPIM     ,JPRB


IMPLICIT NONE

! Declaration of arguments


REAL(KIND=JPRB)             , INTENT(OUT) :: PNORM(:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PSPEC(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSET(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KMASTER
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PMET(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KRESOL

!     ------------------------------------------------------------------

END SUBROUTINE SPECNORM

END INTERFACE
