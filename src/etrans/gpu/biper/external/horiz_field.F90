SUBROUTINE HORIZ_FIELD(KX,KY,PHFIELD)

!   purpose  :
!   --------
!    To produce test horizontal field of temperature.

!    method  :
!   ---------
!    Test horizontal input field is on horizontal grid size KXxKY points, and it
!    represent's temperature. It is obtained form flollwing expression:
!    PHFIELD(i,j)=280*(1+0.1*Sin[PPI*(i+0.5*IMAX)*(j+0.7*IMAX)/IMAX^2+1]) (Pierre Benard) 

!   interface  :
!   ---------
!    CALL HORIZ_FIELD(KX,KY,PHFIELD)
 
!   Explicit arguments :
!   -------------------
!    KX         - number of grid points in x
!    KY         - number of grid points in y
!    PHFIELD    - simulated 2D temperature horizontal field  

!   externals :
!   ----------
!    None.

!   references :
!   ----------

!    author :
!    ------
!    23-May-2008   Antonio Stanesic
!    ----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

!    ----------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),     INTENT(IN)      :: KX
INTEGER(KIND=JPIM),     INTENT(IN)      :: KY
REAL(KIND=JPRB),        INTENT(OUT)     :: PHFIELD(KX,KY)

!    ----------------------------------------------------------------------

REAL(KIND=JPRB),        PARAMETER       :: PPI=3.141592
INTEGER(KIND=JPIM)                      :: JX,JY,IMAX
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!    ----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('HORIZ_FIELD',0,ZHOOK_HANDLE)
!    ----------------------------------------------------------------------

IMAX=MAX(KX,KY)

DO JY=1,KY
 DO JX=1,KX
  PHFIELD(JX,JY)=280*(1+0.1*SIN(PPI*(JX+0.5*IMAX)*(JY+0.7*IMAX)/(IMAX**2)+1))
 ENDDO
ENDDO

!    ----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('HORIZ_FIELD',1,ZHOOK_HANDLE)
END SUBROUTINE HORIZ_FIELD
