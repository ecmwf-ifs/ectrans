INTERFACE
SUBROUTINE EGPNORM_TRANS(PGP,KFIELDS,KPROMA,PAVE,PMIN,PMAX,LDAVE_ONLY,KRESOL)


!**** *EGPNORM_TRANS* - calculate grid-point norms

!     Purpose.
!     --------
!        calculate grid-point norms using a 2 stage (NPRTRV,NPRTRW) communication rather
!        than an approach using a more expensive global gather collective communication

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
!        A.Bogatchev after gpnorm_trans

!     Modifications.
!     --------------
!        Original : 12th Jun 2009

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

! Declaration of arguments
  
REAL(KIND=JPRB),INTENT(IN)    :: PGP(:,:,:)
REAL(KIND=JPRB),INTENT(OUT)   :: PAVE(:)
REAL(KIND=JPRB),INTENT(INOUT) :: PMIN(:)
REAL(KIND=JPRB),INTENT(INOUT) :: PMAX(:)
INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS
INTEGER(KIND=JPIM),INTENT(IN) :: KPROMA
LOGICAL,INTENT(IN)            :: LDAVE_ONLY
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KRESOL

END SUBROUTINE EGPNORM_TRANS
END INTERFACE
