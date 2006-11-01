MODULE SUMPLATBEQ_MOD
CONTAINS
SUBROUTINE SUMPLATBEQ(KDGSA,KDGL,KPROC,KPROCA,KLOENG,LDSPLIT,LDEQ_REGIONS,&
                    &KMEDIAP,KRESTM,KINDIC,KLAST)

!**** *SUMPLATBEQ * - Routine to initialize parallel environment
!                     (latitude partitioning for LEQ_REGIONS=T)

!     Purpose.
!     --------


!**   Interface.
!     ----------
!        *CALL* *SUMPLATBEQ *

!     Explicit arguments - input :
!     -------------------- 
!                          KDGSA      -first latitude (grid-space)
!                                      (may be different from NDGSAG)
!                          KDGL       -last  latitude
!                          KPROC      -total number of processors
!                          KPROCA     -number of processors in A direction
!                          KLOENG     -actual number of longitudes per latitude.
!                          LDSPLIT    -true for latitudes shared between sets
!                          LDEQ_REGIONS -true if eq_regions partitioning

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
!        G. Mozdzynski

!     Modifications.
!     --------------
!        Original : April 2006
!     ------------------------------------------------------------------


USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DISTR
USE EQ_REGIONS_MOD
USE ABORT_TRANS_MOD

IMPLICIT NONE


!     * DUMMY:
INTEGER(KIND=JPIM),INTENT(IN)  :: KDGSA
INTEGER(KIND=JPIM),INTENT(IN)  :: KDGL
INTEGER(KIND=JPIM),INTENT(IN)  :: KPROC
INTEGER(KIND=JPIM),INTENT(IN)  :: KPROCA
INTEGER(KIND=JPIM),INTENT(IN)  :: KLOENG(KDGSA:KDGL)
LOGICAL,INTENT(IN)  :: LDSPLIT
LOGICAL,INTENT(IN)  :: LDEQ_REGIONS
INTEGER(KIND=JPIM),INTENT(OUT)  :: KMEDIAP
INTEGER(KIND=JPIM),INTENT(OUT)  :: KRESTM
INTEGER(KIND=JPIM),INTENT(OUT)  :: KINDIC(KPROCA)
INTEGER(KIND=JPIM),INTENT(OUT)  :: KLAST(KPROCA)

!     * LOCAL:

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: ICOMP, IGL, IMAXI, IMEDIA, IMEDIAP, ITOT, JA, JB, IA, JGL,&
            &ILAST,IREST,IPE,I2REGIONS
LOGICAL   :: LLDONE

!      -----------------------------------------------------------------

!*       1.    COMPUTATION OF KMEDIAP, KRESTM, KINDIC, KLAST.
!              ----------------------------------------------

!     * Computation of KMEDIAP and KRESTM.

IMEDIA = SUM(KLOENG(KDGSA:KDGL))
KMEDIAP = IMEDIA / KPROC

IF( KPROC > 1 )THEN
! test if KMEDIAP is too small and no more than 2 asets would be required 
! for the first latitude
  IF( LDSPLIT )THEN
    I2REGIONS=N_REGIONS(1)+N_REGIONS(2)
    IF( KMEDIAP < (KLOENG(KDGSA)-1)/I2REGIONS+1 )THEN
      WRITE(0,'("SUMPLATBEQ: KMEDIAP=",I6," I2REGIONS=",I3," KLOENG(KDGSA)=",I3)')&
      &KMEDIAP,I2REGIONS,KLOENG(KDGSA)
      CALL ABORT_TRANS ('SUMPLATBEQ: NPROC TOO BIG FOR THIS RESOLUTION, LDSPLIT=T')
    ENDIF
  ELSE
! test for number asets too large for the number of latitudes
    IF( KPROCA > KDGL )THEN
      WRITE(0,'("SUMPLATBEQ: KMEDIAP=",I6," KPROCA=",I4," KDGL=",I4)')&
      &KMEDIAP,KPROCA,KDGL
      CALL ABORT_TRANS ('SUMPLATBEQ: NPROC TOO BIG FOR THIS RESOLUTION, LDSPLIT=F')
    ENDIF
  ENDIF
ENDIF

KRESTM = IMEDIA - KMEDIAP * KPROC
IF (KRESTM  >  0) KMEDIAP = KMEDIAP + 1

!     * Computation of intermediate quantities KINDIC and KLAST

IF (LDSPLIT) THEN

  IREST = 0
  ILAST =0
  IPE=0
  DO JA=1,KPROCA
    ICOMP=0
    DO JB=1,N_REGIONS(JA)
      IPE=IPE+1
      IF (IPE  <=  KRESTM .OR. KRESTM  ==  0) THEN
        ICOMP = ICOMP + KMEDIAP
      ELSE
        ICOMP = ICOMP + (KMEDIAP-1)
      ENDIF
    ENDDO
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
  LLDONE=.FALSE.
  IMEDIAP=KMEDIAP
  IF( MYPROC == 1 )THEN
    WRITE(0,'("SUMPLATBEQ: IMEDIAP=",I6)')IMEDIAP
  ENDIF
  DO WHILE(.NOT.LLDONE)
!   loop until a satisfactory distribution can be found
    IA=1
    IMAXI=IMEDIAP*N_REGIONS(IA)
    DO JGL=1,KDGL
      KLAST(IA)=JGL
      IMAXI=IMAXI-KLOENG(JGL)
      IF( IA == KPROCA .AND. JGL == KDGL )THEN
        IF( MYPROC == 1 )THEN
          WRITE(0,'("SUMPLATBEQ: EXIT 1")')
        ENDIF
        EXIT
      ENDIF
      IF( IA == KPROCA .AND. JGL < KDGL )THEN
        IF( MYPROC == 1 )THEN
          WRITE(0,'("SUMPLATBEQ: EXIT 2")')
        ENDIF
        KLAST(KPROCA)=KDGL
        EXIT
      ENDIF
      IF( IA < KPROCA .AND. JGL == KDGL )THEN
        DO JA=KPROCA,IA+1,-1
          KLAST(JA)=KDGL+JA-KPROCA
        ENDDO
        DO JA=KPROCA,2,-1
          IF( KLAST(JA) <= KLAST(JA-1) )THEN
            KLAST(JA-1)=KLAST(JA)-1
          ENDIF
        ENDDO
        IF( MYPROC == 1 )THEN
          WRITE(0,'("SUMPLATBEQ: EXIT 3")')
        ENDIF
        EXIT
      ENDIF
      IF( IMAXI <= 0 )THEN
        IA=IA+1
        IMAXI=IMAXI+IMEDIAP*N_REGIONS(IA)
      ENDIF
    ENDDO
    IF( KPROCA > 1 .AND. KLAST(KPROCA) == KLAST(KPROCA-1) )THEN
      IMEDIAP=IMEDIAP-1
      IF( MYPROC == 1 )THEN
        WRITE(0,'("SUMPLATBEQ: REDUCING IMEDIAP=",I6)')IMEDIAP
      ENDIF
      IF( IMEDIAP <= 0 )THEN
        CALL ABORT_TRANS ('SUMPLATBEQ: PROBLEM PARTITIONING WITH LSPLIT=F, IMEDIAP <= 0')
      ENDIF
    ELSE
      LLDONE=.TRUE.
    ENDIF
  ENDDO
    
ENDIF

END SUBROUTINE SUMPLATBEQ
END MODULE SUMPLATBEQ_MOD
