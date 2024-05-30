MODULE SUEMPLAT_MOD
CONTAINS
SUBROUTINE SUEMPLAT(KDGL,KPROC,KPROCA,KMYSETA,LDSPLIT,LDEQ_REGIONS,&
 & KFRSTLAT,KLSTLAT,KFRSTLOFF,KPTRLAT,&
 & KPTRFRSTLAT,KPTRLSTLAT,KPTRFLOFF,&
 & PWEIGHT,LDWEIGHTED_DISTR,PMEDIAP,KPROCAGP,&
 & KMEDIAP,KRESTM,LDSPLITLAT,KMYPROC,KLOEN,KDGUX)

!**** *SUEMPLAT * - Initialize gridpoint distrbution in N-S direction

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *SUEMPLAT *

!     Explicit arguments - input :
!     --------------------
!                          KDGL       -last  latitude
!                          KPROC      -total number of processors
!                          KPROCA     -number of processors in A direction
!                          KMYSETA    -process number in A direction
!                          LDSPLIT    -true for latitudes shared between sets
!                          PWEIGHT    -weight per grid-point if weighted
!                                      distribution
!                          LDEQ_REGIONS -true if eq_regions partitioning
!                          LDWEIGHTED_DISTR -true if weighted distribution

!     Explicit arguments - output:
!     --------------------
!                          KMEDIAP    -mean number of grid points per PE
!                          KRESTM     -number of PEs with one extra point
!                          KFRSTLAT   -first latitude row on processor
!                          KLSTLAT    -last  latitude row on processor
!                          KFRSTLOFF  -offset for first latitude in set
!                          KPROCAGP   -number of grid points per A set
!                          KPTRLAT    -pointer to start of latitude
!                          KPTRFRSTLAT-pointer to first latitude
!                          KPTRLSTLAT -pointer to last  latitude
!                          KPTRFLOFF  -offset for pointer to first latitude
!                          LDSPLITLAT -true for latitudes which are split
!                          PMEDIAP    -mean weight per PE if weighted
!                                          distribution
!

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.   SUMPLATB and SUEMPLATB.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        MPP Group *ECMWF*

!     Modifications.
!     --------------
!        Original : 95-10-01
!        David Dent:97-06-02 parameters KFRSTLAT etc added
!        JF. Estrade:97-11-13 Adaptation to ALADIN case
!        J.Boutahar: 98-07-06  phasing with CY19
!        Modified 98-08-10 by K. YESSAD: removal of LRPOLE option + cleanings
!         (correct computation of extrapolar latitudes for KPROCL).
!        Modified 98-12-07 by K. YESSAD and C. FISCHER: cleaning.
!         - merge old sumplat.F and suemplat.F
!         - gather 'lelam' code and 'not lelam' code.
!         - clean (useless duplication of variables, non doctor features).
!         - remodularise according to lelam/not lelam
!           -> lelam features in new routine suemplatb.F,
!              not lelam features in new routine sumplatb.F
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        A.Bogatchev   20-Sep-2010 Phasing cy37
!        R. El Khatib 09-Aug-2013 Allow LEQ_REGIONS
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_GEN         ,ONLY : NOUT, NPRINTLEV

USE SUEMPLATB_MOD   ,ONLY : SUEMPLATB
USE SUMPLATBEQ_MOD  ,ONLY : SUMPLATBEQ
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
!

IMPLICIT NONE

!     * DUMMY:
INTEGER(KIND=JPIM),INTENT(OUT) :: KMEDIAP
INTEGER(KIND=JPIM),INTENT(OUT) :: KRESTM
INTEGER(KIND=JPIM),INTENT(IN)  :: KDGL
INTEGER(KIND=JPIM),INTENT(IN)  :: KPROC
INTEGER(KIND=JPIM),INTENT(IN)  :: KPROCA
INTEGER(KIND=JPIM),INTENT(IN)  :: KMYSETA
INTEGER(KIND=JPIM),INTENT(OUT) :: KFRSTLAT(:)
INTEGER(KIND=JPIM),INTENT(OUT) :: KLSTLAT(:)
INTEGER(KIND=JPIM),INTENT(OUT) :: KFRSTLOFF
INTEGER(KIND=JPIM),INTENT(OUT) :: KPTRLAT(:)
INTEGER(KIND=JPIM),INTENT(OUT) :: KPTRFRSTLAT(:)
INTEGER(KIND=JPIM),INTENT(OUT) :: KPTRLSTLAT(:)
INTEGER(KIND=JPIM),INTENT(OUT) :: KPTRFLOFF
INTEGER(KIND=JPIM),INTENT(OUT) :: KPROCAGP(KPROCA)
REAL(KIND=JPRB),INTENT(OUT)    :: PMEDIAP
REAL(KIND=JPRB),INTENT(IN)     :: PWEIGHT(:)

LOGICAL,INTENT(IN)  :: LDSPLIT
LOGICAL,INTENT(IN)  :: LDEQ_REGIONS
LOGICAL,INTENT(OUT) :: LDSPLITLAT(:)
LOGICAL,INTENT(INOUT)          :: LDWEIGHTED_DISTR
INTEGER(KIND=JPIM),INTENT(IN) :: KMYPROC
INTEGER(KIND=JPIM),INTENT(IN) :: KLOEN(KDGL)
INTEGER(KIND=JPIM),INTENT(IN) :: KDGUX

! === END OF INTERFACE BLOCK ===
INTEGER(KIND=JPIM) :: INDIC(KPROCA),ILAST(KPROCA)

INTEGER(KIND=JPIM) ::  IPTRLATITUDE, JA, JGL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

LOGICAL :: LLFOURIER
LOGICAL  :: LLDEBUG=.FALSE.

!      -----------------------------------------------------------------

!*       1.    CODE DEPENDING ON 'LELAM': COMPUTATION OF
!              KMEDIAP, KRESTM, INDIC, ILAST.
!              -----------------------------------------

IF (LHOOK) CALL DR_HOOK('SUEMPLAT_MOD:SUEMPLAT',0,ZHOOK_HANDLE)

INDIC(:)=0
ILAST(:)=0

IF(LDWEIGHTED_DISTR.AND..NOT.LDEQ_REGIONS)THEN
  CALL ABORT_TRANS ('SUEMPLAT: LDWEIGHTED_DISTR=T AND  LDEQ_REGIONS=F NOT SUPPORTED')
ENDIF

IF( LDEQ_REGIONS )THEN
  CALL SUMPLATBEQ(1,KDGL,KPROC,KPROCA,KLOEN,LDSPLIT,LDEQ_REGIONS,&
   &PWEIGHT,LDWEIGHTED_DISTR,PMEDIAP,KPROCAGP,&
   &KMEDIAP,KRESTM,INDIC,ILAST)
ELSE
  LLFOURIER=.FALSE.
!REK  commented out for now ... monkey business to be done again, should lead to the use of sumplatb
!REK  CALL SUMPLATB(1,KDGL,KPROCA,G%NLOEN,LDSPLIT,LLFOURIER,&
!REK   &KMEDIAP,KRESTM,INDIC,ILAST)
  CALL SUEMPLATB(1,KDGL,KPROCA,KLOEN,LDSPLIT,&
   & PWEIGHT,LDWEIGHTED_DISTR,PMEDIAP,KPROCAGP,&
   & KMEDIAP,KRESTM,INDIC,ILAST,KDGUX)
ENDIF

!      -----------------------------------------------------------------

!*       2.    CODE NOT DEPENDING ON 'LELAM': COMPUTATION OF
!              KFRSTLAT TO LDSPLITLAT.
!              ---------------------------------------------

!     * Computation of first and last latitude of processor sets
!       -----------  in grid-point-space -----------------------
IF(KMYPROC==1.AND.LLDEBUG)THEN
  WRITE(0,'("")')
  WRITE(0,'("SUEMPLAT_MOD:LDWEIGHTED_DISTR=",L1)')LDWEIGHTED_DISTR
  WRITE(0,'("")')
  DO JA=1,KPROCA
    WRITE(0,'("SUEMPLAT_MOD: JA=",I3," ILAST=",I3," INDIC=",I3)')&
     &JA,ILAST(JA),INDIC(JA)
  ENDDO
  WRITE(0,'("")')
  IF( LDEQ_REGIONS .AND. LDSPLIT )THEN
    DO JA=1,KPROCA
      WRITE(0,'("SUEMPLAT_MOD: JA=",I3," KPROCAGP=",I8)')&
       &JA,KPROCAGP(JA)
    ENDDO
    WRITE(0,'("")')
  ENDIF
ENDIF
KFRSTLAT(1) = 1
KLSTLAT(KPROCA) = KDGL
DO JA=1,KPROCA-1
  IF(KMYPROC==1 .AND. NPRINTLEV > 1)THEN
    WRITE(NOUT,'("SUEMPLAT_MOD: JA=",I3," ILAST=",I3," INDIC=",I3)')&
    &JA,ILAST(JA),INDIC(JA)
  ENDIF
  IF ((.NOT. LDSPLIT) .OR. INDIC(JA) == 0) THEN
    KFRSTLAT(JA+1) = ILAST(JA) + 1
    KLSTLAT(JA) = ILAST(JA)
  ELSE
    KFRSTLAT(JA+1) = INDIC(JA)
    KLSTLAT(JA) = INDIC(JA)
  ENDIF
ENDDO
KFRSTLOFF=KFRSTLAT(KMYSETA)-1

!     * Initialise following data structures:-
!       NPTRLAT     (pointer to the start of each latitude)
!       LSPLITLAT   (TRUE if latitude is split over two A sets)
!       NPTRFRSTLAT (pointer to the first latitude of each A set)
!       NPTRLSTLAT  (pointer to the last  latitude of each A set)

DO JGL=1,KDGL
  KPTRLAT  (JGL)=-999
  LDSPLITLAT(JGL)=.FALSE.
ENDDO
IPTRLATITUDE=0
DO JA=1,KPROCA
  DO JGL=KFRSTLAT(JA),KLSTLAT(JA)
    IPTRLATITUDE=IPTRLATITUDE+1
    LDSPLITLAT(JGL)=.TRUE.
    IF( KPTRLAT(JGL) == -999 )THEN
      KPTRLAT(JGL)=IPTRLATITUDE
      LDSPLITLAT(JGL)=.FALSE.
    ENDIF
  ENDDO
ENDDO
DO JA=1,KPROCA
  IF( LDSPLITLAT(KFRSTLAT(JA)) .AND. JA /= 1 )THEN
    KPTRFRSTLAT(JA)=KPTRLAT(KFRSTLAT(JA))+1
  ELSE
    KPTRFRSTLAT(JA)=KPTRLAT(KFRSTLAT(JA))
  ENDIF
  IF( LDSPLITLAT(KLSTLAT(JA)) .AND. JA == KPROCA)THEN
    KPTRLSTLAT(JA)=KPTRLAT(KLSTLAT(JA))+1
  ELSE
    KPTRLSTLAT(JA)=KPTRLAT(KLSTLAT(JA))
  ENDIF
ENDDO
KPTRFLOFF=KPTRFRSTLAT(KMYSETA)-1
IF(KMYPROC==1 .AND. NPRINTLEV > 1)THEN
  DO JGL=1,KDGL       
    WRITE(NOUT,'("SUEMPLAT_MOD: JGL=",I3," KPTRLAT=",I3," LDSPLITLAT=",L4)')&
    & JGL,KPTRLAT(JGL),LDSPLITLAT(JGL)
  ENDDO
  DO JA=1,KPROCA
    WRITE(NOUT,'("SUEMPLAT_MOD: JA=",I3," KFRSTLAT=",I3," KLSTLAT=",I3,&
    & " KPTRFRSTLAT=",I3," KPTRLSTLAT=",I3)')&
    & JA,KFRSTLAT(JA),KLSTLAT(JA),KPTRFRSTLAT(JA),KPTRLSTLAT(JA)
  ENDDO
ENDIF

IF (LHOOK) CALL DR_HOOK('SUEMPLAT_MOD:SUEMPLAT',1,ZHOOK_HANDLE)
END SUBROUTINE SUEMPLAT
END MODULE SUEMPLAT_MOD

