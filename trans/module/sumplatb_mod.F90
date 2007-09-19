MODULE SUMPLATB_MOD
CONTAINS
SUBROUTINE SUMPLATB(KDGSA,KDGL,KPROCA,KLOENG,LDSPLIT,LDFOURIER,&
                    &KMEDIAP,KRESTM,KINDIC,KLAST)

!**** *SUMPLATB * - Routine to initialize parallel environment

!     Purpose.
!     --------


!**   Interface.
!     ----------
!        *CALL* *SUMPLATB *

!     Explicit arguments - input :
!     -------------------- 
!                          KDGSA      -first latitude (grid-space)
!                                      (may be different from NDGSAG)
!                          KDGL       -last  latitude
!                          KPROCA     -number of processors in A direction
!                          KLOENG     -actual number of longitudes per latitude.
!                          LDSPLIT    -true for latitudes shared between sets
!                          LDFOURIER  -true for fourier space partitioning

!     Explicit arguments - output:
!     -------------------- 
!                          KMEDIAP    -mean number of grid points per PE
!                          KRESTM     -number of PEs with one extra point
!                          KINDIC     -intermediate quantity for 'sumplat'
!                          KLAST      -intermediate quantity for 'sumplat'

!        Implicit arguments :
!        --------------------


!     Method.
!     -------
!        See documentation

!     Externals.   NONE.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        K. YESSAD (after old version of sumplat.F).

!     Modifications.
!     --------------
!        Original : 98-12-07
!     ------------------------------------------------------------------


USE PARKIND1  ,ONLY : JPIM, JPIB, JPRB

USE ABORT_TRANS_MOD

IMPLICIT NONE


!     * DUMMY:
INTEGER(KIND=JPIM),INTENT(IN)  :: KDGSA
INTEGER(KIND=JPIM),INTENT(IN)  :: KDGL
INTEGER(KIND=JPIM),INTENT(IN)  :: KPROCA
INTEGER(KIND=JPIM),INTENT(IN)  :: KLOENG(KDGSA:KDGL)
LOGICAL,INTENT(IN)  :: LDSPLIT
LOGICAL,INTENT(IN)  :: LDFOURIER
INTEGER(KIND=JPIM),INTENT(OUT)  :: KMEDIAP
INTEGER(KIND=JPIM),INTENT(OUT)  :: KRESTM
INTEGER(KIND=JPIM),INTENT(OUT)  :: KINDIC(KPROCA)
INTEGER(KIND=JPIM),INTENT(OUT)  :: KLAST(KPROCA)

!     * LOCAL:
INTEGER(KIND=JPIM) :: IPP1(KPROCA),ILAST1(KPROCA)
INTEGER(KIND=JPIM) :: IPP(KPROCA)
INTEGER(KIND=JPIM) :: IFIRST(KPROCA)
INTEGER(KIND=JPIB) :: ICOST(KDGSA:KDGL)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: ICOMP, IGL, IMAXI, IMAXIOL, JA, JGL,&
            &ILAST,IREST,ILIMIT,IFRST
INTEGER(KIND=JPIB) :: IMEDIA,ITOT
REAL(KIND=JPRB) :: ZLG
LOGICAL   :: LLDONE

!      -----------------------------------------------------------------

!*       1.    COMPUTATION OF KMEDIAP, KRESTM, KINDIC, KLAST.
!              ----------------------------------------------

!     * Computation of KMEDIAP and KRESTM.

IF( LDFOURIER )THEN

  DO JGL=1,KDGL
    ZLG=LOG(FLOAT(KLOENG(JGL)))
    ICOST(JGL)=KLOENG(JGL)*ZLG*SQRT(ZLG)
  ENDDO

ELSE

  DO JGL=1,KDGL
    ICOST(JGL)=KLOENG(JGL)
  ENDDO

ENDIF
  
IMEDIA = SUM(ICOST(KDGSA:KDGL))
KMEDIAP = IMEDIA / KPROCA
KRESTM = IMEDIA - KMEDIAP * KPROCA
IF (KRESTM  >  0) KMEDIAP = KMEDIAP + 1

!     * Computation of intermediate quantities KINDIC and KLAST

IF (LDSPLIT) THEN

  IREST = 0
  ILAST =0
  DO JA=1,KPROCA
    IF (JA  <=  KRESTM .OR. KRESTM  ==  0) THEN
      ICOMP = KMEDIAP
    ELSE
      ICOMP = KMEDIAP - 1
    ENDIF
    ITOT = IREST
    IGL = ILAST+1
    DO JGL=IGL,KDGL
      ILAST = JGL
      IF(ITOT+ICOST(JGL) < ICOMP) THEN
        ITOT = ITOT+ICOST(JGL)
      ELSEIF(ITOT+ICOST(JGL) == ICOMP) THEN
        IREST = 0
        KLAST(JA) = JGL 
        KINDIC(JA) = 0
        EXIT
      ELSE
        IREST =  ICOST(JGL) -(ICOMP-ITOT)
        KLAST(JA) = JGL 
        KINDIC(JA) = JGL
        EXIT
      ENDIF
    ENDDO
  ENDDO
  
ELSE

  KINDIC(:) = 0

  IMAXI = KMEDIAP-1
  IMAXIOL = HUGE(IMAXIOL)
  DO
    ILIMIT = IMAXI
    IMAXI = 0
    IFRST = KDGL
    ILAST1(:) = 0
    IPP1(:) = 0
    DO JA=KPROCA,1,-1
      IGL = IFRST
      LATS:DO JGL=IGL,1,-1
        IF (IPP1(JA) < ILIMIT .OR. JA == 1) THEN
          IFRST = JGL-1
          IPP1(JA) = IPP1(JA) + ICOST(JGL)
          IF(ILAST1(JA)  ==  0) ILAST1(JA) = JGL
        ELSE
          EXIT LATS
        ENDIF
      ENDDO LATS
      IMAXI = MAX (IMAXI,IPP1(JA))
    ENDDO
    IF(IMAXI >= IMAXIOL) EXIT
    KLAST(:) = ILAST1(:)
    IPP(:) = IPP1(:)
    IMAXIOL = IMAXI
  ENDDO
  
!       make the distribution more uniform
!       ----------------------------------
  
  IFIRST(1) = 0
  IF (KLAST(1) > 0) IFIRST(1) = 1
  DO JA=2,KPROCA
    IF (IPP(JA) > 0) THEN
      IFIRST(JA) = KLAST(JA-1)+1
    ELSE
      IFIRST(JA) = 0
    ENDIF
  ENDDO
  
  LLDONE = .FALSE.
  DO WHILE( .NOT.LLDONE )
    LLDONE = .TRUE.
  
    DO JA=1,KPROCA-1
      IF (IPP(JA) > IPP(JA+1)) THEN
        IF (IPP(JA)-IPP(JA+1)  >  IPP(JA+1) + 2 *&
         &ICOST(KLAST(JA)) -IPP(JA) ) THEN
          IPP(JA) = IPP(JA) - ICOST(KLAST(JA))
          IPP(JA+1) = IPP(JA+1) + ICOST(KLAST(JA))
          IF (KLAST(JA+1)  ==  0) KLAST(JA+1) = KLAST(JA)
          IFIRST(JA+1) = KLAST(JA)
          KLAST(JA) = KLAST(JA) - 1
          IF (KLAST(JA) == 0) IFIRST(JA) = 0
          LLDONE = .FALSE.
        ENDIF
      ELSE
        IF( IFIRST(JA+1) > 0 )THEN
          IF (IPP(JA+1)-IPP(JA)  >=  IPP(JA) + 2 *&
           &ICOST(IFIRST(JA+1)) -IPP(JA+1) ) THEN
            IPP(JA) = IPP(JA) + ICOST(IFIRST(JA+1))
            IPP(JA+1) = IPP(JA+1) - ICOST(IFIRST(JA+1))
            KLAST(JA) = IFIRST(JA+1)
            IF (IFIRST(JA) == 0) IFIRST(JA) = KLAST(JA)
            IF (KLAST(JA+1)  ==  KLAST(JA)) THEN
              KLAST(JA+1) = 0
              IFIRST(JA+1) = 0
            ELSE
              IFIRST(JA+1) = IFIRST(JA+1) + 1
            ENDIF
           LLDONE = .FALSE.
          ENDIF
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  
ENDIF
  
END SUBROUTINE SUMPLATB
END MODULE SUMPLATB_MOD
