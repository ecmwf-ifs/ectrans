module sumplat_mod
contains
SUBROUTINE SUMPLAT(KDGL,KPROCA,KMYSETA,LDSPLIT,&
                   &KFRSTLAT,KLSTLAT,KFRSTLOFF,KPTRLAT,&
                   &KPTRFRSTLAT,KPTRLSTLAT,KPTRFLOFF,&
                   &KMEDIAP,KRESTM,LDSPLITLAT)

!**** *SUMPLAT * - Initialize gridpoint distrbution in N-S direction

!     Purpose.
!     --------


!**   Interface.
!     ----------
!        *CALL* *SUMPLAT *

!     Explicit arguments - input :
!     -------------------- 
!                          KDGL       -last  latitude
!                          KPROCA     -number of processors in A direction
!                          KMYSETA    -process number in A direction
!                          LDSPLIT    -true for latitudes shared between sets

!     Explicit arguments - output:
!     -------------------- 
!                          KMEDIAP    -mean number of grid points per PE
!                          KRESTM     -number of PEs with one extra point
!                          KFRSTLAT   -first latitude row on processor
!                          KLSTLAT    -last  latitude row on processor
!                          KFRSTLOFF  -offset for first latitude in set
!                          KPTRLAT    -pointer to start of latitude    
!                          KPTRFRSTLAT-pointer to first latitude
!                          KPTRLSTLAT -pointer to last  latitude
!                          KPTRFLOFF  -offset for pointer to first latitude
!                          LDSPLITLAT -true for latitudes which are split

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
!     ------------------------------------------------------------------

#include "tsmbkind.h"

use tpm_geometry

use sumplatb_mod

IMPLICIT NONE


!     * DUMMY:
INTEGER_M,INTENT(OUT) :: KMEDIAP
INTEGER_M,INTENT(OUT) :: KRESTM
INTEGER_M,INTENT(IN)  :: KDGL
INTEGER_M,INTENT(IN)  :: KPROCA
INTEGER_M,INTENT(IN)  :: KMYSETA
INTEGER_M,INTENT(OUT) :: KFRSTLAT(:)
INTEGER_M,INTENT(OUT) :: KLSTLAT(:)
INTEGER_M,INTENT(OUT) :: KFRSTLOFF
INTEGER_M,INTENT(OUT) :: KPTRLAT(:)
INTEGER_M,INTENT(OUT) :: KPTRFRSTLAT(:)
INTEGER_M,INTENT(OUT) :: KPTRLSTLAT(:)
INTEGER_M,INTENT(OUT) :: KPTRFLOFF
LOGICAL,INTENT(IN)  :: LDSPLIT
LOGICAL,INTENT(OUT) :: LDSPLITLAT(:)

!     * LOCAL:
! === END OF INTERFACE BLOCK ===
INTEGER_M :: INDIC(KPROCA),ILAST(KPROCA)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IA, ILAT, IPTRLATITUDE, ISTART, JA, JGL, JGLSUR, JLTLOC


!      -----------------------------------------------------------------

!*       1.    CODE DEPENDING ON 'LELAM': COMPUTATION OF
!              KMEDIAP, KRESTM, INDIC, ILAST.
!              -----------------------------------------


CALL SUMPLATB(1,KDGL,KPROCA,G%NLOEN,LDSPLIT,&
 &KMEDIAP,KRESTM,INDIC,ILAST)

!      -----------------------------------------------------------------

!*       2.    CODE NOT DEPENDING ON 'LELAM': COMPUTATION OF
!              KFRSTLAT TO LDSPLITLAT.
!              ---------------------------------------------


!     * Computation of first and last latitude of processor sets
!       -----------  in grid-point-space -----------------------

KFRSTLAT(1) = 1
KLSTLAT(KPROCA) = KDGL
DO JA=1,KPROCA-1
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
  IF( LDSPLITLAT(KFRSTLAT(JA)) )THEN
    KPTRFRSTLAT(JA)=KPTRLAT(KFRSTLAT(JA))+1
  ELSE
    KPTRFRSTLAT(JA)=KPTRLAT(KFRSTLAT(JA))
  ENDIF
  KPTRLSTLAT(JA)=KPTRLAT(KLSTLAT(JA))
ENDDO
KPTRFLOFF=KPTRFRSTLAT(KMYSETA)-1

RETURN
END SUBROUTINE SUMPLAT
end module sumplat_mod
