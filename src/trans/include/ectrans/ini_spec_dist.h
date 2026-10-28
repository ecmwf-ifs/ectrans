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
!> This subroutine initializes the arrays controlling spectral wave distribution.
!!
!! @param KSMAX Spectral truncation required.
!! @param KTMAX Overtruncation for KSMAX.
!! @param KPRTRW Number of processors in A-direction.
!! @param KMYSETW A-set for present processor.
!! @param KASM0 Offsets for spectral waves.
!! @param KSPOLEGL Local version of NSPOLEG.
!! @param KPROCM Where a certain spectral wave belongs.
!! @param KUMPP Number of spectral waves on this task.
!! @param KSPEC Local version of NSPEC.
!! @param KSPEC2 Local version of NSPEC2.
!! @param KSPEC2MX Maximum KSPEC2 across tasks.
!! @param KPOSSP Global spectral fields partitioning.
!! @param KMYMS This task's spectral zonal wavenumbers.
!! @param KPTRMS Pointer to the first wave number of a given A-set.
!! @param KALLMS Wave numbers for all wave-set concatenated together to give all wave numbers in
!!   wave-set order.
SUBROUTINE INI_SPEC_DIST(KSMAX,KTMAX,KPRTRW,KMYSETW,KASM0,KSPOLEGL,KPROCM,&
                    &KUMPP,KSPEC,KSPEC2,KSPEC2MX,KPOSSP,KMYMS,KPTRMS,KALLMS)


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
