MODULE SUEMPLATB_MOD
CONTAINS
SUBROUTINE SUEMPLATB(KDGSA,KDGL,KPROCA,KLOENG,LDSPLIT,&
 & PWEIGHT,LDWEIGHTED_DISTR,PMEDIAP,KPROCAGP,&
 & KMEDIAP,KRESTM,KINDIC,KLAST,KDGUX)  

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
!                          KDGUX      -last latitude for meaningful computations
!                                      (suggested to pass NDGUX in gp-space, NDGL in Fourier space
!                                       for having a good load-balance)
!                          PWEIGHT    -weight per grid-point if weighted distribution
!                          LDWEIGHTED_DISTR -true if weighted distribution`

!     Explicit arguments - output:
!     -------------------- 
!                          KMEDIAP    -mean number of grid points per PE
!                          KPROCAGP   -number of grid points per A set
!                          KRESTM     -number of PEs with one extra point
!                          KINDIC     -intermediate quantity for 'sumplat'
!                          KLAST      -intermediate quantity for 'sumplat'
!                          PMEDIAP    -mean weight per PE if weighted distribution

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
!         G. Radnoti: 03-03-03: Semi-merge with sumplatb, only difference:
!                               NS-partitioning according to NDGUX
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        A.Bogatchev   21-Sep-2010 phasing CY37
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
!

IMPLICIT NONE

!     * DUMMY:
INTEGER(KIND=JPIM),INTENT(IN)  :: KDGSA
INTEGER(KIND=JPIM),INTENT(IN)  :: KDGL
INTEGER(KIND=JPIM),INTENT(IN)  :: KPROCA
INTEGER(KIND=JPIM),INTENT(IN)  :: KLOENG(KDGSA:KDGL)
REAL(KIND=JPRB),INTENT(IN)     :: PWEIGHT(:)
LOGICAL,INTENT(IN)  :: LDSPLIT
LOGICAL,INTENT(INOUT)  :: LDWEIGHTED_DISTR
INTEGER(KIND=JPIM),INTENT(IN)  :: KDGUX
INTEGER(KIND=JPIM),INTENT(OUT)  :: KMEDIAP
INTEGER(KIND=JPIM),INTENT(OUT)  :: KRESTM
INTEGER(KIND=JPIM),INTENT(OUT)  :: KINDIC(KPROCA)
INTEGER(KIND=JPIM),INTENT(OUT)  :: KLAST(KPROCA)
INTEGER(KIND=JPIM),INTENT(IN)  :: KPROCAGP(KPROCA)
REAL(KIND=JPRB),INTENT(IN)     :: PMEDIAP

INTEGER(KIND=JPIM) :: IPP1(KPROCA),ILAST1(KPROCA)
INTEGER(KIND=JPIM) :: IPP(KPROCA)
INTEGER(KIND=JPIM) :: IFIRST(KPROCA)

INTEGER(KIND=JPIM) :: ICOMP, IGL, IMAXI, IMAXIOL, IMEDIA, ITOT, JA, JGL,&
 & ILAST,IREST,ILIMIT,IFRST
LOGICAL   :: LLDONE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      -----------------------------------------------------------------

!*       1.    COMPUTATION OF KMEDIAP, KRESTM, KINDIC, KLAST.
!              ----------------------------------------------

!     * Computation of KMEDIAP and KRESTM.

IF (LHOOK) CALL DR_HOOK('SUEMPLATB_MOD:SUEMPLATB',0,ZHOOK_HANDLE)
IF (LDWEIGHTED_DISTR) THEN
  CALL ABORT_TRANS ('SUMPLATBEQ: ALADIN CODE IS NOT PREPARED FOR WEIGHTED DISTRIBUTION')
ENDIF
IMEDIA = SUM(KLOENG(KDGSA:KDGUX))
KMEDIAP = IMEDIA / KPROCA
IF (KMEDIAP  <  KLOENG(KDGL/2)) THEN
  CALL ABORT_TRANS ('SUMPLATB: KPROCA TOO BIG FOR THIS RESOLUTION')
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
    DO JGL=IGL,KDGUX
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
  KLAST(KPROCA)=KDGL
  KINDIC(KPROCA)=0
ELSE

  KINDIC(:) = 0

  IMAXI = KMEDIAP-1
  IMAXIOL = HUGE(IMAXIOL)
  DO
    ILIMIT = IMAXI
    IMAXI = 0
    IFRST = KDGUX
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
           & KLOENG(KLAST(JA)) -IPP(JA) ) THEN  
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
             & KLOENG(IFIRST(JA+1)) -IPP(JA+1) ) THEN  
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
  KLAST(KPROCA)=KDGL
ENDIF

IF (LHOOK) CALL DR_HOOK('SUEMPLATB_MOD:SUEMPLATB',1,ZHOOK_HANDLE)
END SUBROUTINE SUEMPLATB
END MODULE SUEMPLATB_MOD
