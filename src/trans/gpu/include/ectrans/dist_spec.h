INTERFACE
SUBROUTINE DIST_SPEC(PSPECG,KFDISTG,KFROM,KVSET,KRESOL,PSPEC,&
 & LDIM1_IS_FLD,KSMAX,KSORT)

!**** *DIST_SPEC* - Distribute global spectral array among processors

!     Purpose.
!     --------
!        Interface routine for distributing spectral array

!**   Interface.
!     ----------
!     CALL DIST__SPEC(...)

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

!     Externals.  SET_RESOL   - set resolution
!     ----------  DIST_SPEC_CONTROL - control routine

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
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPEC(:,:)
LOGICAL            ,OPTIONAL, INTENT(IN)  :: LDIM1_IS_FLD
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KSMAX
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KSORT (:)


!     ------------------------------------------------------------------

END SUBROUTINE DIST_SPEC

END INTERFACE
