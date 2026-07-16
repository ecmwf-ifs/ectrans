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
!                    partition; the grid parameters are then unused.
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
!    PMU(KSIZEJ), PGW(KSIZEJ) : sin(latitude) and Gaussian weights (global)

USE ISO_FORTRAN_ENV, ONLY: INT64
USE PARKIND1, ONLY: JPIM, JPRB, JPRD

IMPLICIT NONE

INTEGER(KIND=INT64), INTENT(OUT) :: KRETURNCODE
INTEGER(KIND=INT64), INTENT(IN)  :: KRESOL, KSIZEJ, KTRUNC, KSLOEN, KNUMMAXRESOL
INTEGER(KIND=INT64), DIMENSION(KSLOEN), INTENT(IN) :: KLOEN
INTEGER(KIND=INT64), INTENT(OUT) :: KGPTOT, KSPEC, KSPEC2, KGPTOTG, KSPEC2G, KSMAX
INTEGER(KIND=INT64), DIMENSION(KSIZEJ), INTENT(OUT) :: KNMENG
REAL(KIND=JPRB), DIMENSION(KSIZEJ), INTENT(OUT) :: PMU, PGW
!
! II. Local variables declaration
INTEGER(KIND=JPIM) :: IRESOL, IIDENTRESOL
INTEGER(KIND=JPIM) :: IGPTOT, ISPEC, ISPEC2, IGPTOTG, ISPEC2G, ISMAX
INTEGER(KIND=JPIM) :: INMENG(KSIZEJ)
INTEGER, DIMENSION(KSLOEN) :: ILOEN
REAL(KIND=JPRD) :: ZMU(KSIZEJ)
REAL(KIND=JPRB) :: ZDELTAX, ZDELTAY
LOGICAL :: LLSTOP
#include "trans_inq.h"
KRETURNCODE = 0
LLSTOP = .FALSE.
ZMU(:) = 0._JPRD
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
  IRESOL = INT(KRESOL, JPIM)
ENDIF
!
IF (.NOT. LLSTOP) THEN
  CALL TRANS_INQ(KRESOL=IRESOL, KGPTOT=IGPTOT, KSPEC=ISPEC, KSPEC2=ISPEC2, &
   & KGPTOTG=IGPTOTG, KSPEC2G=ISPEC2G, KSMAX=ISMAX, KNMENG=INMENG, PMU=ZMU, PGW=PGW)
  KGPTOT  = IGPTOT
  KSPEC   = ISPEC
  KSPEC2  = ISPEC2
  KGPTOTG = IGPTOTG
  KSPEC2G = ISPEC2G
  KSMAX   = ISMAX
  KNMENG(:) = INMENG(:)
  PMU(:) = ZMU(:)
ENDIF
END SUBROUTINE TRANS_INQ4PY
