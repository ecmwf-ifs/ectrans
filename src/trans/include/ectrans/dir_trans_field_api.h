! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 

INTERFACE

SUBROUTINE DIR_TRANS_FIELD_API(YDFSPVOR,YDFSPDIV,YDFSPSCALAR, &
    & YDFU, YDFV, YDFSCALAR, &
    & KSPEC, KPROMA, KGPBLKS, KGPTOT, KFLEVG, KFLEVL,KPROC,&
    & LDACC)


!**** *DIR_TRANS_FIELD_API* - Field API interface to direct spectral transform

!     Purpose.
!     --------
!        Allow to call DIR_TRANS with a list of fields from field API                            

!**   Interface.
!     ----------
!     CALL DIR_TRANS_FIELD_API(...)

!     Explicit arguments :
!     --------------------
!      output
!       YDFSPVOR(:)    - List of spectral vector fields (vorticity) 
!       YDFSPDIV(:)    - List of spectral vector fields (divergence)
!       YDFSPSCALAR(:) - List of spectral scalar fields 
!      input
!       YDFU(:)        - List of grid-point vector fields (u)
!       YDFV(:)        - List of grid-point vector fields (v)
!       YDFSCALAR(:)   - List of grid-point scalar fields
!       KSPEC          - Number of spectral coefficients
!       KPROMA         - Blocking factor
!       KGPBLKS        - Number of blocks
!       KGPTOT         - Number of total grid points
!       KFLEVG         - Number of levels
!       KFLEVL         - Number of local levels
!       KPROC          - Processor ID
!       LDACC          - Field data on device

USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE FIELD_API_BASIC_TYPE_MOD, ONLY: FIELD_BASIC_PTR
USE PARKIND1  ,ONLY : JPIM     ,JPRB

TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFSPVOR(:), YDFSPDIV(:)
TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFSPSCALAR(:)

TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFU(:),YDFV(:)
TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFSCALAR(:)

INTEGER(KIND=JPIM), INTENT(IN) ::KSPEC
INTEGER(KIND=JPIM), INTENT(IN) ::KPROMA
INTEGER(KIND=JPIM), INTENT(IN) ::KGPBLKS
INTEGER(KIND=JPIM), INTENT(IN) ::KGPTOT
INTEGER(KIND=JPIM), INTENT(IN) :: KFLEVG
INTEGER(KIND=JPIM), INTENT(IN) :: KFLEVL
INTEGER(KIND=JPIM), INTENT(IN) :: KPROC
LOGICAL, INTENT(IN), OPTIONAL  :: LDACC


END SUBROUTINE DIR_TRANS_FIELD_API
END INTERFACE
