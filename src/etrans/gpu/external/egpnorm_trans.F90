SUBROUTINE EGPNORM_TRANS(PGP,KFIELDS,KPROMA,PAVE,PMIN,PMAX,LDAVE_ONLY,KRESOL)


!**** *EGPNORM_TRANS* - calculate grid-point norms

!     Purpose.
!     --------
!        calculate grid-point norms

!**   Interface.
!     ----------
!     CALL EGPNORM_TRANS(...)

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
!     PMIN        - minimum (input/output)
!     PMAX        - maximum (input/output)
!     LDAVE_ONLY  - T : PMIN and PMAX already contain local MIN and MAX
!     KRESOL      -  resolution tag (optional)
!                    default assumes first defined resolution
!

!     Author.
!     -------
!        George Mozdzynski *ECMWF*

!     Modifications.
!     --------------
!        Original : 19th Sept 2008
!        R. El Khatib 07-08-2009 Optimisation directive for NEC
!        R. El Khatib 16-Sep-2019 merge with global model code
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD

!ifndef INTERFACE

USE TPM_DIM         ,ONLY : R
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_DIM         ,ONLY : R
USE ESET_RESOL_MOD  ,ONLY : ESET_RESOL
USE YOMHOOK         ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!endif INTERFACE

IMPLICIT NONE

! Declaration of arguments

REAL(KIND=JPRB)   ,INTENT(IN)    :: PGP(:,:,:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAVE(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PMIN(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PMAX(:)
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIELDS
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
LOGICAL           ,INTENT(IN)    :: LDAVE_ONLY
INTEGER(KIND=JPIM),OPTIONAL, INTENT(IN)  :: KRESOL

!ifndef INTERFACE

! Local variables
INTEGER(KIND=JPIM) :: JGL
REAL(KIND=JPRB) :: ZW(R%NDGL)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('EGPNORM_TRANS',0,ZHOOK_HANDLE)

! Set current resolution
CALL ESET_RESOL(KRESOL)

DO JGL=1,R%NDGL
  ZW(1:)=1._JPRB/G%NLOEN(JGL)
ENDDO

CALL GPNORM_TRANS_GPU(PGP,KFIELDS,KPROMA,PAVE,PMIN,PMAX,LDAVE_ONLY)

! note: weighting not taken into account by GPNORM_TRANS_GPU, so we do it here
PAVE=PAVE/G%NLOEN(1)

IF (LHOOK) CALL DR_HOOK('EGPNORM_TRANS',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!endif INTERFACE


END SUBROUTINE EGPNORM_TRANS
