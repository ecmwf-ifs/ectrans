! (C) Copyright 2008- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE GPNORM_TRANSTL(PGP,KFIELDS,KPROMA,PAVE,KRESOL)


!**** *GPNORM_TRANSTL* - calculate grid-point norms
!                        reduced version for linear model

!     Purpose.
!     --------
!        calculate grid-point norms

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
!        Filip Vana, (c) ECMWF
!        9-Sep-2024

!     Modifications.
!     --------------

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

!ifndef INTERFACE

USE TPM_DIM         ,ONLY : R
USE TPM_DISTR       ,ONLY : D
USE TPM_FIELDS      ,ONLY : F
USE TPM_DIM         ,ONLY : R
USE SET_RESOL_MOD   ,ONLY : SET_RESOL
USE YOMHOOK         ,ONLY : LHOOK,   DR_HOOK,   JPHOOK
USE GPNORM_TRANS_CTLTL_MOD, ONLY : GPNORM_TRANS_CTLTL

!endif INTERFACE

IMPLICIT NONE

! Declaration of arguments

REAL(KIND=JPRB)   ,INTENT(IN)    :: PGP(:,:,:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAVE(:)
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIELDS
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
INTEGER(KIND=JPIM),OPTIONAL, INTENT(IN)  :: KRESOL

!ifndef INTERFACE

! Local variables
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GPNORM_TRANSTL',0,ZHOOK_HANDLE)

! Set current resolution
CALL SET_RESOL(KRESOL)

CALL GPNORM_TRANS_CTLTL(PGP,KFIELDS,KPROMA,PAVE,F%RW(1:R%NDGL))

IF (LHOOK) CALL DR_HOOK('GPNORM_TRANSTL',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!endif INTERFACE


END SUBROUTINE GPNORM_TRANSTL
