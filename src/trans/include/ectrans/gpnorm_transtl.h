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
SUBROUTINE GPNORM_TRANSTL(PGP,KFIELDS,KPROMA,PAVE,KRESOL)


!**** *GPNORM_TRANSTL* - calculate grid-point norms

!     Purpose.
!     --------
!        calculate grid-point norms using a 2 stage (NPRTRV,NPRTRW) communication rather
!        than an approach using a more expensive global gather collective communication

!**   Interface.
!     ----------
!     CALL GPNORM_TRANSTL(...)

!     Explicit arguments :
!     --------------------
!     PGP(:,:,:) - gridpoint fields (input)
!                  PGP is  dimensioned (NPROMA,KFIELDS,NGPBLKS) where
!                  NPROMA is the blocking factor, KFIELDS the total number
!                  of fields and NGPBLKS the number of NPROMA blocks.
!     KFIELDS     - number of fields (input)
!                   (these do not have to be just levels)
!     KPROMA      - required blocking factor (input)
!     PAVE        - average (output)
!     KRESOL      -  resolution tag (optional)
!                    default assumes first defined resolution
!

!     Author.
!     -------
!        George Mozdzynski *ECMWF*

!     Modifications.
!     --------------
!        Original : 19th Sept 2008

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

! Declaration of arguments
  
REAL(KIND=JPRB),INTENT(IN)    :: PGP(:,:,:)
REAL(KIND=JPRB),INTENT(OUT)   :: PAVE(:)
INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS
INTEGER(KIND=JPIM),INTENT(IN) :: KPROMA
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KRESOL

END SUBROUTINE GPNORM_TRANSTL
END INTERFACE
