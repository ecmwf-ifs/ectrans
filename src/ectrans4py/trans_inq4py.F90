SUBROUTINE TRANS_INQ4PY(KRETURNCODE, KRESOL, KSIZEJ, KTRUNC, KSLOEN, KLOEN, KNUMMAXRESOL, &
                      & KGPTOT, KSPEC, KSPEC2, KGPTOTG, KSPEC2G, KSMAX, KNMENG, PMU, PGW)
! ** PURPOSE
!    Wrapper to TRANS_INQ: extract the geometry of a resolution -- this task's
!    LOCAL sizes (KGPTOT, KSPEC, KSPEC2), the GLOBAL sizes (KGPTOTG, KSPEC2G), the
!    truncation (KSMAX), the cut-off zonal wavenumber per latitude (KNMENG), and the
!    Gaussian latitudes PMU=sin(latitude) and weights PGW (both global, length KSIZEJ).
!
!    Serves both the serial and the distributed setups through KRESOL:
!    * KRESOL <= 0 : self-initialise serially (LDMPOFF) from the grid parameters
!                    (KSIZEJ, KTRUNC, KSLOEN, KLOEN, KNUMMAXRESOL) via SPEC_SETUP4PY,
!                    then inquire. Local == global (single task).
!    * KRESOL >  0 : inquire the resolution already set up by SETUP_TRANS_4PY (after
!                    the parallel SETUP_TRANS0_4PY). Local sizes are this task's
!                    partition; the grid parameters are then unused. A KRESOL that is
!                    not a resolution set up by SETUP_TRANS_4PY returns error code 101,
!                    a grid-point-only one (LDGRIDONLY), without spectral sizes, error
!                    code 102 (TRANS_INQ_GRID4PY inquires its grid).
!
! ** AUTHOR
!    9 April 2014, S. Riette
!
! ** MODIFICATIONS
!    6 Jan., S. Riette: w_spec_setup interfaced modified
!
! I. Dummy arguments declaration
!    KRESOL       : resolution tag from SETUP_TRANS_4PY, or <=0 to self-initialise
!    KSIZEJ       : number of Gaussian latitudes (sizes KNMENG, PMU, PGW)
!    KTRUNC       : spectral truncation             (self-init only)
!    KSLOEN       : size of KLOEN                    (self-init only)
!    KLOEN        : number of points on each latitude row (self-init only)
!    KNUMMAXRESOL : maximum number of resolutions    (self-init only)
!    KGPTOT/KSPEC/KSPEC2      : this task's local grid / spectral / doubled-spectral sizes
!    KGPTOTG/KSPEC2G          : global grid / doubled-spectral sizes
!    KSMAX                    : spectral truncation T
!    KNMENG(KSIZEJ)           : cut-off zonal wavenumber per latitude
!    PMU(KSIZEJ), PGW(KSIZEJ) : sin(latitude) and Gaussian weights (global); PMU in double
!                               precision whatever JPRB is, as TRANS_INQ returns it

USE ISO_FORTRAN_ENV, ONLY: INT64
USE PARKIND1, ONLY: JPIM, JPRB, JPRD
USE RESOL_CHECK4PY_MOD, ONLY: RESOL_STATUS

IMPLICIT NONE

INTEGER(KIND=INT64), INTENT(OUT) :: KRETURNCODE
INTEGER(KIND=INT64), INTENT(IN)  :: KRESOL, KSIZEJ, KTRUNC, KSLOEN, KNUMMAXRESOL
INTEGER(KIND=INT64), DIMENSION(KSLOEN), INTENT(IN) :: KLOEN
INTEGER(KIND=INT64), INTENT(OUT) :: KGPTOT, KSPEC, KSPEC2, KGPTOTG, KSPEC2G, KSMAX
INTEGER(KIND=INT64), DIMENSION(KSIZEJ), INTENT(OUT) :: KNMENG
REAL(KIND=JPRD), DIMENSION(KSIZEJ), INTENT(OUT) :: PMU
REAL(KIND=JPRB), DIMENSION(KSIZEJ), INTENT(OUT) :: PGW
!
! II. Local variables declaration
INTEGER(KIND=JPIM) :: IRESOL, IIDENTRESOL
INTEGER(KIND=JPIM) :: IGPTOT, ISPEC, ISPEC2, IGPTOTG, ISPEC2G, ISMAX
INTEGER(KIND=JPIM) :: INMENG(KSIZEJ)
INTEGER, DIMENSION(KSLOEN) :: ILOEN
REAL(KIND=JPRB) :: ZDELTAX, ZDELTAY
LOGICAL :: LLSTOP
#include "trans_inq.h"
KRETURNCODE = 0
LLSTOP = .FALSE.
INMENG(:) = 0
!
IF (KRESOL <= 0_INT64) THEN
  ! Serial self-initialisation from the grid parameters (LDMPOFF; single task).
  ILOEN(:) = KLOEN(:)
  ZDELTAX = 0.0_JPRB
  ZDELTAY = 0.0_JPRB
  CALL SPEC_SETUP4PY(KRETURNCODE, 0, INT(KSIZEJ), 0, 0, &
                   & INT(KTRUNC), 0, INT(KNUMMAXRESOL), ILOEN, .FALSE., INT(KSLOEN), &
                   & ZDELTAX, ZDELTAY, IIDENTRESOL, LLSTOP)
  IRESOL = IIDENTRESOL
ELSE
  ! Inquire the resolution already set up by SETUP_TRANS_4PY (distributed).
  KRETURNCODE = RESOL_STATUS(KRESOL)
  IF (KRETURNCODE /= 0) RETURN
  IRESOL = INT(KRESOL, JPIM)
ENDIF
!
IF (.NOT. LLSTOP) THEN
  CALL TRANS_INQ(KRESOL=IRESOL, KGPTOT=IGPTOT, KSPEC=ISPEC, KSPEC2=ISPEC2, &
   & KGPTOTG=IGPTOTG, KSPEC2G=ISPEC2G, KSMAX=ISMAX, KNMENG=INMENG, PMU=PMU, PGW=PGW)
  KGPTOT  = IGPTOT
  KSPEC   = ISPEC
  KSPEC2  = ISPEC2
  KGPTOTG = IGPTOTG
  KSPEC2G = ISPEC2G
  KSMAX   = ISMAX
  KNMENG(:) = INMENG(:)
ENDIF
END SUBROUTINE TRANS_INQ4PY


SUBROUTINE TRANS_INQ_VSET4PY(KRETURNCODE, KRESOL, KPRTRW, KMYSETW, KMYSETV)
! ** PURPOSE
!    Wrapper to TRANS_INQ for the spectral decomposition of a distributed resolution, as
!    needed by the V-set (*_VSET4PY) routines: the number of wave sets KPRTRW, and this
!    task's wave set KMYSETW and V-set KMYSETV (both 1-based). The number of V-sets is
!    NPRTRV = NPROC/KPRTRW. A KRESOL that is not a resolution set up by SETUP_TRANS_4PY
!    returns error code 101, a grid-point-only one (LDGRIDONLY) error code 102.
USE ISO_FORTRAN_ENV, ONLY: INT64
USE PARKIND1, ONLY: JPIM
USE RESOL_CHECK4PY_MOD, ONLY: RESOL_STATUS
IMPLICIT NONE
INTEGER(KIND=INT64), INTENT(OUT) :: KRETURNCODE
INTEGER(KIND=INT64), INTENT(IN)  :: KRESOL
INTEGER(KIND=INT64), INTENT(OUT) :: KPRTRW, KMYSETW, KMYSETV
INTEGER(KIND=JPIM) :: IPRTRW, IMYSETW, IMYSETV
#include "trans_inq.h"
KRETURNCODE = RESOL_STATUS(KRESOL)
IF (KRETURNCODE /= 0) RETURN
CALL TRANS_INQ(KRESOL=INT(KRESOL,JPIM), KPRTRW=IPRTRW, KMYSETW=IMYSETW, KMYSETV=IMYSETV)
KPRTRW  = IPRTRW
KMYSETW = IMYSETW
KMYSETV = IMYSETV
END SUBROUTINE TRANS_INQ_VSET4PY


SUBROUTINE TRANS_INQ_GRID4PY(KRETURNCODE, KRESOL, KDGL, KNRNS, KNREW, &
                           & KGPTOT, KGPTOTG, KGPTOTMX, KFRSTLOFF, KPTRFLOFF, &
                           & KMY_REGION_NS, KMY_REGION_EW, KFRSTLAT, KLSTLAT, KPTRFRSTLAT, &
                           & KSTA, KONL, KGPTOTL, PMU)
! ** PURPOSE
!    Wrapper to TRANS_INQ for the grid-point partition tables of a distributed resolution,
!    spectral or grid-point-only (SETUP_TRANS_4PY LDGRIDONLY), as the IFS inquires them
!    for its radiation grid: this task's region, the latitudes of each N-S band (a-set),
!    and the start and length of every latitude row in every region (NSTA, NONL).
!
! ** ARGUMENTS
!    KRESOL        : resolution tag from SETUP_TRANS_4PY
!    KDGL          : number of latitude rows of that resolution
!    KNRNS, KNREW  : numbers of N-S bands and of E-W regions, KREGIONS_NS and KREGIONS_EW
!                    of SETUP_TRANS0_4PY (or SETUP_TRANS0_REGIONS4PY)
! ** RETURNS
!    KGPTOT, KGPTOTG, KGPTOTMX : grid points on this task, on the globe, on any task at most
!    KFRSTLOFF     : offset of the first latitude of this task's N-S band (NFRSTLOFF)
!    KPTRFLOFF     : offset of this band's first row in NSTA/NONL (NPTRFLOFF)
!    KMY_REGION_NS, KMY_REGION_EW : this task's region (1-based)
!    KFRSTLAT(KNRNS), KLSTLAT(KNRNS) : first and last latitude of each N-S band
!    KPTRFRSTLAT(KNRNS) : row of the first latitude of each N-S band in NSTA/NONL
!    KSTA(KDGL+KNRNS-1, KNREW) : first grid point of each row in each region (NSTA): a
!                                latitude split between two bands has a row in each
!    KONL(KDGL+KNRNS-1, KNREW) : number of grid points of each row in each region (NONL)
!    KGPTOTL(KNRNS, KNREW)     : grid points of each region
!    PMU(KDGL)     : sines of the latitudes (global), in double precision
!
! ** RETURN CODE
!    0 success; 1 KNRNS/KNREW differ from the region counts of SETUP_TRANS0; 2 KDGL
!    differs from the number of latitude rows of KRESOL; 101 KRESOL is not a resolution set
!    up by SETUP_TRANS_4PY. They are checked before the table inquiry, which would abort in
!    TRANS_INQ.
USE ISO_FORTRAN_ENV, ONLY: INT64
USE PARKIND1, ONLY: JPIM, JPRD
USE RESOL_CHECK4PY_MOD, ONLY: RESOL_STATUS, REGION_COUNTS
IMPLICIT NONE
INTEGER(KIND=INT64), INTENT(OUT) :: KRETURNCODE
INTEGER(KIND=INT64), INTENT(IN)  :: KRESOL, KDGL, KNRNS, KNREW
INTEGER(KIND=INT64), INTENT(OUT) :: KGPTOT, KGPTOTG, KGPTOTMX, KFRSTLOFF, KPTRFLOFF
INTEGER(KIND=INT64), INTENT(OUT) :: KMY_REGION_NS, KMY_REGION_EW
INTEGER(KIND=INT64), DIMENSION(KNRNS), INTENT(OUT) :: KFRSTLAT, KLSTLAT, KPTRFRSTLAT
INTEGER(KIND=INT64), DIMENSION(KDGL+KNRNS-1, KNREW), INTENT(OUT) :: KSTA, KONL
INTEGER(KIND=INT64), DIMENSION(KNRNS, KNREW), INTENT(OUT) :: KGPTOTL
REAL(KIND=JPRD), DIMENSION(KDGL), INTENT(OUT) :: PMU
INTEGER(KIND=JPIM) :: IRESOL
INTEGER(KIND=INT64) :: IREGIONS_NS, IREGIONS_EW
INTEGER(KIND=JPIM) :: IGPTOT, IGPTOTG, IGPTOTMX, IFRSTLOFF, IPTRFLOFF, IMYNS, IMYEW
INTEGER(KIND=JPIM) :: IFRSTLAT(KNRNS), ILSTLAT(KNRNS), IPTRFRSTLAT(KNRNS)
INTEGER(KIND=JPIM) :: ISTA(KDGL+KNRNS-1, KNREW), IONL(KDGL+KNRNS-1, KNREW)
INTEGER(KIND=JPIM) :: IGPTOTL(KNRNS, KNREW)
#include "trans_inq.h"
KGPTOT = 0
KGPTOTG = 0
KGPTOTMX = 0
KFRSTLOFF = 0
KPTRFLOFF = 0
KMY_REGION_NS = 0
KMY_REGION_EW = 0
KFRSTLAT(:) = 0
KLSTLAT(:) = 0
KPTRFRSTLAT(:) = 0
KSTA(:,:) = 0
KONL(:,:) = 0
KGPTOTL(:,:) = 0
PMU(:) = 0._JPRD
KRETURNCODE = RESOL_STATUS(KRESOL, LDGRIDPOINT_ONLY=.TRUE.)
IF (KRETURNCODE /= 0) RETURN
CALL REGION_COUNTS(IREGIONS_NS, IREGIONS_EW)
IF (KNRNS /= IREGIONS_NS .OR. KNREW /= IREGIONS_EW) THEN
  KRETURNCODE = 1
  RETURN
ENDIF
IRESOL = INT(KRESOL, JPIM)
! The last latitude of the last N-S band is the number of latitude rows
IFRSTLAT(:) = 0
ILSTLAT(:) = 0
CALL TRANS_INQ(KRESOL=IRESOL, KFRSTLAT=IFRSTLAT, KLSTLAT=ILSTLAT)
IF (KDGL /= INT(MAXVAL(ILSTLAT), INT64)) THEN
  KRETURNCODE = 2
  RETURN
ENDIF
ISTA(:,:) = 0
IONL(:,:) = 0
IGPTOTL(:,:) = 0
CALL TRANS_INQ(KRESOL=IRESOL, KGPTOT=IGPTOT, KGPTOTG=IGPTOTG, KGPTOTMX=IGPTOTMX, &
 & KGPTOTL=IGPTOTL, KFRSTLAT=IFRSTLAT, KLSTLAT=ILSTLAT, KFRSTLOFF=IFRSTLOFF, &
 & KPTRFRSTLAT=IPTRFRSTLAT, KPTRFLOFF=IPTRFLOFF, KSTA=ISTA, KONL=IONL, &
 & KMY_REGION_NS=IMYNS, KMY_REGION_EW=IMYEW, PMU=PMU)
KGPTOT = IGPTOT
KGPTOTG = IGPTOTG
KGPTOTMX = IGPTOTMX
KFRSTLOFF = IFRSTLOFF
KPTRFLOFF = IPTRFLOFF
KMY_REGION_NS = IMYNS
KMY_REGION_EW = IMYEW
KFRSTLAT(:) = IFRSTLAT(:)
KLSTLAT(:) = ILSTLAT(:)
KPTRFRSTLAT(:) = IPTRFRSTLAT(:)
KSTA(:,:) = ISTA(:,:)
KONL(:,:) = IONL(:,:)
KGPTOTL(:,:) = IGPTOTL(:,:)
END SUBROUTINE TRANS_INQ_GRID4PY
