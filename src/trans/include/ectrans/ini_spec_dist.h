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
SUBROUTINE INI_SPEC_DIST(KSMAX,KTMAX,KPRTRW,KMYSETW,KASM0,KSPOLEGL,KPROCM,&
                    &KUMPP,KSPEC,KSPEC2,KSPEC2MX,KPOSSP,KMYMS,KPTRMS,KALLMS)


! begin_doc_block
! ## `INI_SPEC_DIST`
!
! ### Signature
!
! ```f90
! SUBROUTINE INI_SPEC_DIST(KSMAX, KTMAX, KPRTRW, KMYSETW, KASM0, KSPOLEGL, KPROCM, KUMPP, KSPEC, &
!   &                      KSPEC2, KSPEC2MX, KPOSSP, KMYMS, KPTRMS, KALLMS)
! ```
!
! ### Purpose
!
! This subroutine initializes the arrays controlling spectral wave distribution.
!
! ### `INTENT(IN)` arguments
!
! - `INTEGER(KIND=JPIM), INTENT(IN) :: KSMAX`  
!   Spectral truncation required.
! - `INTEGER(KIND=JPIM), INTENT(IN) :: KTMAX`  
!   Overtruncation for KSMAX.
! - `INTEGER(KIND=JPIM), INTENT(IN) :: KPRTRW`  
!   Number of processors in A-direction.
! - `INTEGER(KIND=JPIM), INTENT(IN) :: KMYSETW`  
!   A-set for present processor.
!
! ### `OPTIONAL, INTENT(OUT)` arguments
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KASM0`  
!   Offsets for spectral waves.
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KSPOLEGL`  
!   Local version of NSPOLEG.
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KPROCM`  
!   Where a certain spectral wave belongs.
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KUMPP`  
!   Number of spectral waves on this task.
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KSPEC`  
!   Local version of NSPEC.
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KSPEC2`
!   Local version of NSPEC2.
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KSPEC2MX`  
!   Maximum KSPEC2 across tasks.
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KPOSSP`  
!   Global spectral fields partitioning.
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KMYMS`  
!   This task's spectral zonal wavenumbers.
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KPTRMS`  
!   Pointer to the first wave number of a given A-set.
! - `INTEGER(KIND=JPIM), OPTIONAL, INTENT(OUT) :: KALLMS`  
!   Wave numbers for all wave-set concatenated together to give all wave numbers in wave-set order.
! end_doc_block

USE EC_PARKIND, ONLY: JPIM
IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KSMAX
INTEGER(KIND=JPIM),INTENT(IN)  :: KTMAX
INTEGER(KIND=JPIM),INTENT(IN)  :: KPRTRW
INTEGER(KIND=JPIM),INTENT(IN)  :: KMYSETW
INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KSPEC
INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KSPEC2
INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KSPEC2MX
INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KSPOLEGL

INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KASM0(0:KSMAX)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KPROCM(0:KSMAX)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KUMPP(KPRTRW)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KPOSSP(KPRTRW+1)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KMYMS(KSMAX+1)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KPTRMS(KPRTRW)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KALLMS(KSMAX+1)
END SUBROUTINE INI_SPEC_DIST
END INTERFACE
