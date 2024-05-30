INTERFACE
SUBROUTINE EDIST_SPEC(PSPECG,KFDISTG,KFROM,KVSET,KRESOL,PSPEC,&
  & LDIM1_IS_FLD,KSORT)

!**** *EDIST_SPEC* - Distribute global spectral array among processors

!     Purpose.
!     --------
!        Interface routine for distributing spectral array

!**   Interface.
!     ----------
!     CALL EDIST__SPEC(...)

!     Explicit arguments : 
!     -------------------- 
!     PSPECG(:,:) - Global spectral array
!     KFDISTG     - Global number of fields to be distributed
!     KFROM(:)    - Processor resposible for distributing each field
!     KVSET(:)    - "B-Set" for each field
!     KRESOL      - resolution tag  which is required ,default is the
!                   first defined resulution (input)
!     PSPEC(:,:)  - Local spectral array
!
!     Method.
!     -------

!     Externals.  ESET_RESOL   - set resolution
!     ----------  EDIST_SPEC_CONTROL - control routine

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

REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PSPECG(:,:)
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFDISTG
INTEGER(KIND=JPIM)          , INTENT(IN)  :: KFROM(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSET(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KRESOL
LOGICAL            ,OPTIONAL, INTENT(IN)  :: LDIM1_IS_FLD
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPEC(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KSORT (:)


!     ------------------------------------------------------------------

END SUBROUTINE EDIST_SPEC
END INTERFACE
