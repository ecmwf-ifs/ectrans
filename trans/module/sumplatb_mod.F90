module sumplatb_mod
contains
SUBROUTINE SUMPLATB(KDGSA,KDGL,KPROCA,KLOENG,LDSPLIT,&
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


#include "tsmbkind.h"

IMPLICIT NONE


!     * DUMMY:
INTEGER_M,INTENT(IN)  :: KDGSA
INTEGER_M,INTENT(IN)  :: KDGL
INTEGER_M,INTENT(IN)  :: KPROCA
INTEGER_M,INTENT(IN)  :: KLOENG(KDGSA:KDGL)
LOGICAL,INTENT(IN)  :: LDSPLIT
INTEGER_M,INTENT(OUT)  :: KMEDIAP
INTEGER_M,INTENT(OUT)  :: KRESTM
INTEGER_M,INTENT(OUT)  :: KINDIC(KPROCA)
INTEGER_M,INTENT(OUT)  :: KLAST(KPROCA)

!     * LOCAL:
INTEGER_M :: IPP1(KPROCA),ILAST1(KPROCA)
INTEGER_M :: IPP(KPROCA)
INTEGER_M :: IFIRST(KPROCA)

!     LOCAL INTEGER SCALARS
INTEGER_M :: ICOMP, IGL, IMAXI, IMAXIOL, IMEDIA, IND, ITOT, JA, JGL,&
            &ILAST,IREST,ILIMIT,IFRST
LOGICAL   :: LLDONE

!      -----------------------------------------------------------------

!*       1.    COMPUTATION OF KMEDIAP, KRESTM, KINDIC, KLAST.
!              ----------------------------------------------

!     * Computation of KMEDIAP and KRESTM.

IMEDIA = SUM(KLOENG(KDGSA:KDGL))
KMEDIAP = IMEDIA / KPROCA
IF (KMEDIAP  <  KLOENG(KDGL/2)) THEN
  CALL ABOR1 ('SUMPLATB: KPROCA TOO BIG FOR THIS RESOLUTION')
ENDIF
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
      IF(ITOT+KLOENG(JGL) < ICOMP) THEN
        ITOT = ITOT+KLOENG(JGL)
      ELSEIF(ITOT+KLOENG(JGL) == ICOMP) THEN
        IREST = 0
        KLAST(JA) = JGL 
        KINDIC(JA) = 0
        EXIT
      ELSE
        IREST =  KLOENG(JGL) -(ICOMP-ITOT)
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
          IPP1(JA) = IPP1(JA) + KLOENG(JGL)
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
         &KLOENG(KLAST(JA)) -IPP(JA) ) THEN
          IPP(JA) = IPP(JA) - KLOENG(KLAST(JA))
          IPP(JA+1) = IPP(JA+1) + KLOENG(KLAST(JA))
          IF (KLAST(JA+1)  ==  0) KLAST(JA+1) = KLAST(JA)
          IFIRST(JA+1) = KLAST(JA)
          KLAST(JA) = KLAST(JA) - 1
          IF (KLAST(JA) == 0) IFIRST(JA) = 0
          LLDONE = .FALSE.
        ENDIF
      ELSE
        IF( IFIRST(JA+1) > 0 )THEN
          IF (IPP(JA+1)-IPP(JA)  >=  IPP(JA) + 2 *&
           &KLOENG(IFIRST(JA+1)) -IPP(JA+1) ) THEN
            IPP(JA) = IPP(JA) + KLOENG(IFIRST(JA+1))
            IPP(JA+1) = IPP(JA+1) - KLOENG(IFIRST(JA+1))
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

RETURN
END SUBROUTINE SUMPLATB
end module sumplatb_mod
