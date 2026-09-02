MODULE ESEAM_MOD
CONTAINS
SUBROUTINE ESEAM(KDLUN,KDLON,KDGUN,KDGL,KDLUX,KDGUX,KNUBI,PWORK,KSEAM,PSEAM,LDBIX,LDBIY)

!   purpose  :
!   --------
!     Make ESEAM extension.

!    *CALL* *ESEAM*(...)

!      externals :
!      ----------
!             None

!      explicit arguments :
!      ------------------
!     KDLUN : lower bound for the x (or longitude) dimension
!             of the gridpoint array
!     KDLON : upper bound for the x (or longitude) dimension
!             of the gridpoint array on C U I U E
!     KDGUN : lower bound for the y (or latitude) dimension
!             of the gridpoint array
!     KDGL : upper bound for the y (or latitude) dimension
!             of the gridpoint array on C U I U E
!     KDLUX : upper bound for the x (or longitude) dimension
!             of  C U I.
!     KDGUX : upper bound for the y (or latitude) dimension
!             of  C U I.
!     KNUBI  : number of levels to biperiodicise
!     PWORK : gridpoint array on C U I U E.
!     KSEAM : size of the SEAM mixing zone
!     PSEAM : scalar parameter for SEAM (Lmix and Lboyd)
!     LDBIX : .TRUE. biperiodicisation in x  ( and vice versa )
!     LDBIY : .TRUE. biperiodicisation in y  ( and vice versa )

!      references :
!      ----------

!      author :
!      ------
!       B. Menetrier

!      modifications :
!      -------------
!      B. Menetrier : 02-06-2025: Initial implementation
!      -------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK,  JPHOOK

!      -------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KDLUN
INTEGER(KIND=JPIM),INTENT(IN)    :: KDLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGUN
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGL
INTEGER(KIND=JPIM),INTENT(IN)    :: KDLUX
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGUX
INTEGER(KIND=JPIM),INTENT(IN)    :: KNUBI
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWORK(KDLUN:KDLON,KNUBI,KDGUN:KDGL)
INTEGER(KIND=JPIM),INTENT(IN)    :: KSEAM
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSEAM(2)
LOGICAL           ,INTENT(IN)    :: LDBIX
LOGICAL           ,INTENT(IN)    :: LDBIY

!      -------------------------------------------------------------

INTEGER(KIND=JPIM) :: JFL, JLON, JGL, JJ
INTEGER(KIND=JPIM) :: NEXTX, NEXTY, JEXTX, JEXTY
INTEGER(KIND=JPIM) :: IDLEFT, IDRIGHT, IDBOT, IDTOP
INTEGER(KIND=JPIM) :: IWLEFT, IWRIGHT, IWBOT, IWTOP
REAL(KIND=JPRB) :: ZX(2), ZY(2), ZSYM, ZASYM, ZMIX, ZBOYD
REAL(KIND=JPRB), ALLOCATABLE :: ZFACTX(:,:), ZFACTY(:,:)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      -------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ESEAM',0,ZHOOK_HANDLE)
!      -------------------------------------------------------------

! Check if there is anything to do
IF((.NOT.LDBIX).AND.(.NOT.LDBIY)) THEN
  IF (LHOOK) CALL DR_HOOK('ESEAM',1,ZHOOK_HANDLE)
  RETURN
ENDIF

! Extension zone size
NEXTX = KDLON-KDLUX
NEXTY = KDGL-KDGUX

! Allocation
ALLOCATE(ZFACTX(NEXTX,2))
ALLOCATE(ZFACTY(NEXTY,2))

! Mixing / Boyd factors
DO JEXTX=1,NEXTX
  ZX(1) = REAL(JEXTX,JPRB)/REAL(MIN(NEXTX,KSEAM)+1)
  ZX(2) = REAL(JEXTX,JPRB)/REAL(NEXTX+1)
  DO JJ=1,2
    IF (ZX(JJ)<1.0_JPRB) THEN
      ZFACTX(JEXTX,JJ) = 0.5_JPRB*(1.0_JPRB+ERF(PSEAM(JJ)*(1.0_JPRB-2.0_JPRB*ZX(JJ))/SQRT(4.0_JPRB*ZX(JJ)*(1.0_JPRB-ZX(JJ)))))
    ELSE
      ZFACTX(JEXTX,JJ) = 0.0_JPRB
    ENDIF
  ENDDO
ENDDO
DO JEXTY=1,NEXTY
  ZY(1) = REAL(JEXTY,JPRB)/REAL(MIN(NEXTY,KSEAM)+1)
  ZY(2) = REAL(JEXTY,JPRB)/REAL(NEXTY+1)
  DO JJ=1,2
    IF (ZY(JJ)<1.0_JPRB) THEN
      ZFACTY(JEXTY,JJ) = 0.5_JPRB*(1.0_JPRB+ERF(PSEAM(JJ)*(1.0_JPRB-2.0_JPRB*ZY(JJ))/SQRT(4.0_JPRB*ZY(JJ)*(1.0_JPRB-ZY(JJ)))))
    ELSE
      ZFACTY(JEXTY,JJ) = 0.0_JPRB
    ENDIF
  ENDDO
ENDDO

! Loop over levels
DO JFL=1,KNUBI
  DO JGL=KDGUX+1,KDGL
    DO JLON=KDLUX+1,KDLON
      ! Initialize
      PWORK(JLON,JFL,JGL) = 0.0_JPRB

      ! Distances to boundaries
      IDLEFT = NEXTX-(JLON-KDLUX)+1
      IDRIGHT = JLON-KDLUX
      IDBOT = NEXTY-(JGL-KDGUX)+1
      IDTOP = JGL-KDGUX
      IF ((IDLEFT<1).OR.(IDLEFT>NEXTX)) CALL ABOR1("WRONG IDLEFT INDEX")
      IF ((IDRIGHT<1).OR.(IDRIGHT>NEXTX)) CALL ABOR1("WRONG IDRIGHT INDEX")
      IF ((IDBOT<1).OR.(IDBOT>NEXTY)) CALL ABOR1("WRONG IDBOT INDEX")
      IF ((IDTOP<1).OR.(IDTOP>NEXTY)) CALL ABOR1("WRONG IDTOP INDEX")


      ! Work array indices
      IWLEFT = 1+IDLEFT
      IWRIGHT = KDLUX-IDRIGHT
      IWBOT = 1+IDBOT
      IWTOP = KDGUX-IDTOP

      ! Left-bottom corner
      ZSYM = PWORK(IWLEFT,JFL,IWBOT)
      ZASYM = 2.0_JPRB*PWORK(1,JFL,1)-PWORK(IWLEFT,JFL,IWBOT)
      ZMIX = ZFACTX(IDLEFT,1)*ZFACTY(IDBOT,1)
      ZBOYD = ZFACTX(IDLEFT,2)*ZFACTY(IDBOT,2)
      PWORK(JLON,JFL,JGL) = PWORK(JLON,JFL,JGL)+(ZASYM*ZMIX+ZSYM*(1.0_JPRB-ZMIX))*ZBOYD

      ! Left-top corner
      ZSYM = PWORK(IWLEFT,JFL,IWTOP)
      ZASYM = 2.0_JPRB*PWORK(1,JFL,KDGUX)-PWORK(IWLEFT,JFL,IWTOP)
      ZMIX = ZFACTX(IDLEFT,1)*ZFACTY(IDTOP,1)
      ZBOYD = ZFACTX(IDLEFT,2)*ZFACTY(IDTOP,2)
      PWORK(JLON,JFL,JGL) = PWORK(JLON,JFL,JGL)+(ZASYM*ZMIX+ZSYM*(1.0_JPRB-ZMIX))*ZBOYD

      ! Right-top corner
      ZSYM = PWORK(IWRIGHT,JFL,IWTOP)
      ZASYM = 2.0_JPRB*PWORK(KDLUX,JFL,KDGUX)-PWORK(IWRIGHT,JFL,IWTOP)
      ZMIX = ZFACTX(IDRIGHT,1)*ZFACTY(IDTOP,1)
      ZBOYD = ZFACTX(IDRIGHT,2)*ZFACTY(IDTOP,2)
      PWORK(JLON,JFL,JGL) = PWORK(JLON,JFL,JGL)+(ZASYM*ZMIX+ZSYM*(1.0_JPRB-ZMIX))*ZBOYD

      ! Right-bottom corner
      ZSYM = PWORK(IWRIGHT,JFL,IWBOT)
      ZASYM = 2.0_JPRB*PWORK(KDLUX,JFL,1)-PWORK(IWRIGHT,JFL,IWBOT)
      ZMIX = ZFACTX(IDRIGHT,1)*ZFACTY(IDBOT,1)
      ZBOYD = ZFACTX(IDRIGHT,2)*ZFACTY(IDBOT,2)
      PWORK(JLON,JFL,JGL) = PWORK(JLON,JFL,JGL)+(ZASYM*ZMIX+ZSYM*(1.0_JPRB-ZMIX))*ZBOYD
    ENDDO
  ENDDO

  DO JGL=1,KDGL
    DO JLON=KDLUX+1,KDLON
      ! Initialize
      PWORK(JLON,JFL,JGL) = 0.0_JPRB

      ! Distances to boundaries
      IDLEFT = NEXTX-(JLON-KDLUX)+1
      IDRIGHT = JLON-KDLUX
      IF ((IDLEFT<1).OR.(IDLEFT>NEXTX)) CALL ABOR1("WRONG IDLEFT INDEX")
      IF ((IDRIGHT<1).OR.(IDRIGHT>NEXTX)) CALL ABOR1("WRONG IDRIGHT INDEX")

      ! Work array indices
      IWLEFT = 1+IDLEFT
      IWRIGHT = KDLUX-IDRIGHT

      ! Left side
      ZSYM = PWORK(IWLEFT,JFL,JGL)
      ZASYM = 2.0_JPRB*PWORK(1,JFL,JGL)-PWORK(IWLEFT,JFL,JGL)
      ZMIX = ZFACTX(IDLEFT,1)
      ZBOYD = ZFACTX(IDLEFT,2)
      PWORK(JLON,JFL,JGL) = PWORK(JLON,JFL,JGL)+(ZASYM*ZMIX+ZSYM*(1.0_JPRB-ZMIX))*ZBOYD

      ! Right side
      ZSYM = PWORK(IWRIGHT,JFL,JGL)
      ZASYM = 2.0_JPRB*PWORK(KDLUX,JFL,JGL)-PWORK(IWRIGHT,JFL,JGL)
      ZMIX = ZFACTX(IDRIGHT,1)
      ZBOYD = ZFACTX(IDRIGHT,2)
      PWORK(JLON,JFL,JGL) = PWORK(JLON,JFL,JGL)+(ZASYM*ZMIX+ZSYM*(1.0_JPRB-ZMIX))*ZBOYD
    ENDDO
  ENDDO

  DO JGL=KDGUX+1,KDGL
    DO JLON=1,KDLON
      ! Initialize
      PWORK(JLON,JFL,JGL) = 0.0_JPRB

      ! Distances to boundaries
      IDBOT = NEXTY-(JGL-KDGUX)+1
      IDTOP = JGL-KDGUX
      IF ((IDBOT<1).OR.(IDBOT>NEXTY)) CALL ABOR1("WRONG IDBOT INDEX")
      IF ((IDTOP<1).OR.(IDTOP>NEXTY)) CALL ABOR1("WRONG IDTOP INDEX")

      ! Work array indices
      IWBOT = 1+IDBOT
      IWTOP = KDGUX-IDTOP

      ! Bottom side
      ZSYM = PWORK(JLON,JFL,IWBOT)
      ZASYM = 2.0_JPRB*PWORK(JLON,JFL,1)-PWORK(JLON,JFL,IWBOT)
      ZMIX = ZFACTY(IDBOT,1)
      ZBOYD = ZFACTY(IDBOT,2)
      PWORK(JLON,JFL,JGL) = PWORK(JLON,JFL,JGL)+(ZASYM*ZMIX+ZSYM*(1.0_JPRB-ZMIX))*ZBOYD

      ! Top side
      ZSYM = PWORK(JLON,JFL,IWTOP)
      ZASYM = 2.0_JPRB*PWORK(JLON,JFL,KDGUX)-PWORK(JLON,JFL,IWTOP)
      ZMIX = ZFACTY(IDTOP,1)
      ZBOYD = ZFACTY(IDTOP,2)
      PWORK(JLON,JFL,JGL) = PWORK(JLON,JFL,JGL)+(ZASYM*ZMIX+ZSYM*(1.0_JPRB-ZMIX))*ZBOYD
    ENDDO
  ENDDO
ENDDO

! Release memory
DEALLOCATE(ZFACTX)
DEALLOCATE(ZFACTY)

!      -------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ESEAM',1,ZHOOK_HANDLE)
END SUBROUTINE ESEAM
END MODULE ESEAM_MOD
