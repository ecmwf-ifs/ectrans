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
SUBROUTINE EDIST_GRID(PGPG,KPROMA,KFDISTG,KFROM,KRESOL,PGP,KSORT)

!**** *EDIST_GRID* - Distribute global gridpoint array among processors

!     Purpose.
!     --------
!        Interface routine for distributing gridpoint array

!**   Interface.
!     ----------
!     CALL EDIST_GRID(...)

!     Explicit arguments : 
!     -------------------- 
!     PGPG(:,:) - Global spectral array
!     KFDISTG     - Global number of fields to be distributed
!     KPROMA      - required blocking factor for gridpoint input
!     KFROM(:)    - Processor resposible for distributing each field
!     KRESOL      - resolution tag  which is required ,default is the
!                   first defined resulution (input)
!     PGP(:,:)  - Local spectral array
!
!     Method.
!     -------

!     Externals.  ESET_RESOL      - set resolution
!     ----------  DIST_GRID_CTL  - control routine

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03

!     ------------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

! Declaration of arguments

REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGPG(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KPROMA
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFDISTG
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFROM(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KRESOL
REAL(KIND=JPRB)             , INTENT(OUT) :: PGP(:,:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KSORT (:)


!     ------------------------------------------------------------------

END SUBROUTINE EDIST_GRID
END INTERFACE
