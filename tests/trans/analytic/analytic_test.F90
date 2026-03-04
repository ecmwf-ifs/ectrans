! (C) Copyright 2023- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

PROGRAM ANALYTIC_TEST

USE PARKIND1, ONLY: JPIM, JPRB, JPRD, JPRM
USE MPL_MODULE, ONLY: MPL_INIT, MPL_NPROC, MPL_MYRANK, MPL_END
USE ANALYTIC_UTILS, ONLY: ANALYTIC_INIT, PREPARE_LEGENDRE_POLYNOMIALS, COMPUTE_ANALYTIC_SOLUTION, COMPUTE_MAX_ERROR, &
  &                       ANALYTIC_END

IMPLICIT NONE

! Output unit numbers
INTEGER(KIND=JPIM), PARAMETER :: NERR     = 0 ! Unit number for stderr
INTEGER(KIND=JPIM), PARAMETER :: NOUT     = 6 ! Unit number for stdout

! Default parameters
INTEGER(KIND=JPIM) :: NSMAX   = 21  ! Spectral truncation
INTEGER(KIND=JPIM) :: NFLD    = 1   ! Number of scalar fields
LOGICAL            :: LIMAG = .FALSE. ! Test imaginary part of spectral data
INTEGER(KIND=JPIM) :: NDGL ! Number of latitudes
INTEGER(KIND=JPIM) :: NSPEC2
INTEGER(KIND=JPIM) :: NGPTOT
INTEGER(KIND=JPIM) :: NGPTOTG
INTEGER(KIND=JPIM) :: NSPEC2G
INTEGER(KIND=JPIM) :: JA
INTEGER(KIND=JPIM) :: IB
INTEGER(KIND=JPIM) :: JPRTRV
INTEGER(KIND=JPIM) :: N_REGIONS_NS
INTEGER(KIND=JPIM) :: N_REGIONS_EW

INTEGER(KIND=JPIM), ALLOCATABLE :: NLOEN(:)
INTEGER(KIND=JPIM) :: MYPROC

! Spectral and grid point arrays
REAL(KIND=JPRB), ALLOCATABLE :: ZSPSCALAR(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZGP(:,:,:)

! Spectral space data structures
REAL(KIND=JPRD), ALLOCATABLE :: ZSPH_ANALYTIC(:,:)

LOGICAL :: LUSEFLT = .FALSE. ! Use fast legendre transforms
REAL(KIND=JPRD) :: RTOLERANCE = 1E-9 ! maximum relative lmax error tolerance for
                                     ! passing analytyic solution tests in double precision
                                     ! Value for single precision is set at the beginning
                                     ! of the execution.

! Verbosity level (-1, 0 or 1)
INTEGER :: VERBOSITY = -1

INTEGER(KIND=JPIM) :: NPROC ! Number of procs
INTEGER(KIND=JPIM) :: NPRGPNS ! Grid-point decomp
INTEGER(KIND=JPIM) :: NPRGPEW ! Grid-point decomp
INTEGER(KIND=JPIM) :: NPRTRV = 0 ! Spectral decomp
INTEGER(KIND=JPIM) :: NPRTRW = 0 ! Spectral decomp
INTEGER(KIND=JPIM) :: MYSETV

INTEGER(KIND=JPIM), ALLOCATABLE :: NUM_LOCAL_LEVELS_ALL(:), IVSET(:)

INTEGER(KIND=JPIM) :: NFLD_LOCAL, NUM_REST
INTEGER(KIND=JPIM) :: I
INTEGER(KIND=JPIM) :: ISQR
INTEGER(KIND=JPIM) :: NPROMA = 16
INTEGER(KIND=JPIM) :: NGPBLKS
INTEGER(KIND=JPIM) :: ILEV, JLEV
INTEGER(KIND=JPIM) :: JSETV
INTEGER(KIND=JPIM) :: M, N, IMAG_IDX
INTEGER, EXTERNAL :: EC_MPIRANK
LOGICAL :: LUSE_MPI = .TRUE.
CHARACTER(LEN=16) :: CGRID = ''

REAL(KIND=JPRB) :: ZMAX_ERROR

!===================================================================================================

#include "setup_trans0.h"
#include "setup_trans.h"
#include "inv_trans.h"
#include "dir_trans.h"
#include "trans_inq.h"
#include "trans_end.h"
#include "abor1.intfb.h"

!===================================================================================================

LUSE_MPI = DETECT_MPIRUN()
IF (JPRB == JPRM) RTOLERANCE = 1E-3 ! Tolerance for single precision
! Setup
CALL GET_COMMAND_LINE_ARGUMENTS(NSMAX, CGRID, NFLD, LUSEFLT, NPROMA, VERBOSITY, NPRTRV, NPRTRW, &
  &                             LIMAG, RTOLERANCE)
IF (CGRID == '') CGRID = CUBIC_FULL_GRID(NSMAX)
CALL PARSE_GRID(CGRID, NDGL, NLOEN)

WRITE(NOUT,'(A)') '======= ecTrans analytic test ======='
WRITE(NOUT,'(A,I0)') 'Spectral truncation: ', NSMAX
WRITE(NOUT,'(A,A)') 'Grid: ', TRIM(CGRID)
WRITE(NOUT,'(A,I0)') 'Number of scalar fields: ', NFLD

!===================================================================================================

IF (LUSE_MPI) THEN
  CALL MPL_INIT(LDINFO=(VERBOSITY >= 1))
  NPROC  = MPL_NPROC()
  MYPROC = MPL_MYRANK()
ELSE
  NPROC = 1
  MYPROC = 1
ENDIF

!===================================================================================================

! only output to stdout on pe 1
IF (MYPROC /= 1) THEN
  OPEN(UNIT=NOUT, FILE='/dev/null')
ENDIF

!===================================================================================================

! Compute nprgpns and nprgpew
! This version selects most square-like distribution
ISQR = INT(SQRT(REAL(NPROC, JPRB)))
DO JA = ISQR, NPROC
  IB = NPROC / JA
  IF (JA * IB == NPROC) THEN
    NPRGPNS = MAX(JA, IB)
    NPRGPEW = MIN(JA, IB)
    EXIT
  ENDIF
ENDDO

! Compute nprtrv and nprtrw if not provided on the command line
IF (NPRTRV > 0 .OR. NPRTRW > 0) THEN
  IF (NPRTRV == 0) NPRTRV = NPROC / NPRTRW
  IF (NPRTRW == 0) NPRTRW = NPROC / NPRTRV
  IF (NPRTRW * NPRTRV /= NPROC) CALL ABOR1('ANALYTIC_TEST: NPRTRW * NPRTRV /= NPROC')
ELSE
  DO JPRTRV = 4, NPROC
    NPRTRV = JPRTRV
    NPRTRW = NPROC / NPRTRV
    IF (NPRTRV * NPRTRW /= NPROC) CYCLE
    IF (NPRTRV > NPRTRW) EXIT
  ENDDO
  ! Go for approx square partition for backup
  IF (NPRTRV * NPRTRW /= NPROC .OR. NPRTRV > NPRTRW) THEN
    ISQR = INT(SQRT(REAL(NPROC, JPRB)))
    DO JA = ISQR, NPROC
      IB = NPROC / JA
      IF (JA * IB == NPROC) THEN
        NPRTRW = MAX(JA, IB)
        NPRTRV = MIN(JA, IB)
        EXIT
      ENDIF
    ENDDO
  ENDIF
ENDIF

!===================================================================================================
! Call ecTrans setup routines
!===================================================================================================

CALL SETUP_TRANS0(KOUT=NOUT, KERR=NERR, KPRINTLEV=MERGE(2, 0, VERBOSITY == 1), KPRGPNS=NPRGPNS, &
  &               KPRGPEW=NPRGPEW, KPRTRW=NPRTRW, LDALLOPERM=.TRUE., LDMPOFF=.NOT.LUSE_MPI, &
  &               K_REGIONS_NS=N_REGIONS_NS, K_REGIONS_EW=N_REGIONS_EW)

CALL SETUP_TRANS(KSMAX=NSMAX, KDGL=NDGL, KLOEN=NLOEN, LDUSEFLT=LUSEFLT)

CALL TRANS_INQ(KSPEC2=NSPEC2, KSPEC2G=NSPEC2G, KGPTOT=NGPTOT, KGPTOTG=NGPTOTG)

IF (NPROMA == 0) THEN ! no blocking (default when not specified)
  NPROMA = NGPTOT
ENDIF

! Calculate number of NPROMA blocks
NGPBLKS = (NGPTOT - 1) / NPROMA + 1

IF (NPROC == 1) THEN
  MYSETV = 1
ELSE
  CALL TRANS_INQ(KMYSETV=MYSETV)
ENDIF

! Determine number of local levels for fourier and legendre calculations
! based on the values of nfld and nprtrv
ALLOCATE(NUM_LOCAL_LEVELS_ALL(NPRTRV))
NFLD_LOCAL = NFLD / NPRTRV ! INTEGER DIVISION
NUM_REST = NFLD - NFLD_LOCAL * NPRTRV ! TAIL BLOCK
DO JSETV = 1, NPRTRV
  NUM_LOCAL_LEVELS_ALL(JSETV) = NFLD_LOCAL
  IF (JSETV == NPRTRV) THEN
    NUM_LOCAL_LEVELS_ALL(JSETV) = NUM_LOCAL_LEVELS_ALL(JSETV) + NUM_REST
  ENDIF
ENDDO

NFLD_LOCAL = NUM_LOCAL_LEVELS_ALL(MYSETV)

!===================================================================================================
! Allocate arrays
!===================================================================================================

ALLOCATE(IVSET(NFLD))

! Compute spectral distribution
ILEV = 0
DO JSETV = 1, NPRTRV
  DO JLEV = 1, NUM_LOCAL_LEVELS_ALL(JSETV)
    ILEV = ILEV + 1
    IVSET(ILEV) = JSETV
  ENDDO
ENDDO

ALLOCATE(ZSPSCALAR(NFLD_LOCAL, NSPEC2))
ALLOCATE(ZGP(NPROMA, NFLD, NGPBLKS))

! Allocate arrays for analytic solutions
ALLOCATE(ZSPH_ANALYTIC(NPROMA, NGPBLKS))

! Compute geographic longitude GELAM and latitude GELAT:
CALL ANALYTIC_INIT(NPROMA, NGPBLKS, NDGL, N_REGIONS_NS, N_REGIONS_EW, NLOEN, NSMAX)

!===================================================================================================
! Perform tests
!===================================================================================================

ZMAX_ERROR = 0.0_JPRB

! Loop over all wavenumbers (check actually tested wavenumber inside)
DO M = 0, NSMAX
  ! Prepare Legendre polynomials for this zonal wavenumber
  CALL PREPARE_LEGENDRE_POLYNOMIALS(M, NSMAX)

  DO N = M, NSMAX
    CALL INITIALIZE_SPECTRAL_ARRAY(NSMAX, M, N, ZSPSCALAR)

    ! Compute analytic solutions
    ZSPH_ANALYTIC = COMPUTE_ANALYTIC_SOLUTION(NPROMA, NGPBLKS, NGPTOT, M, N)

    ! Do inverse transform
    CALL INV_TRANS(PSPSCALAR=ZSPSCALAR, KPROMA=NPROMA, KVSETSC=IVSET, PGP=ZGP)

    ZMAX_ERROR = MAX(ZMAX_ERROR, COMPUTE_MAX_ERROR(NGPTOT, NPROMA, ZGP(:, 1, :), ZSPH_ANALYTIC(:, :)))

    ! Do direct transform
    CALL DIR_TRANS(PGP=ZGP, KPROMA=NPROMA, KVSETSC=IVSET, PSPSCALAR=ZSPSCALAR)
  END DO
END DO

!===================================================================================================
! Cleanup
!===================================================================================================

DEALLOCATE(ZSPH_ANALYTIC, ZGP, ZSPSCALAR)

!===================================================================================================
! Finalize
!===================================================================================================

CALL TRANS_END

IF (LUSE_MPI) THEN
  CALL MPL_END(LDMEMINFO=.FALSE.)
ENDIF

!===================================================================================================
! Close file
!===================================================================================================

IF (NPROC > 1) THEN
  IF (MYPROC /= 1) THEN
    CLOSE(UNIT=NOUT)
  ENDIF
ENDIF

IF (ZMAX_ERROR > RTOLERANCE) THEN
  WRITE(NERR, '(A)') '*******************************'
  WRITE(NERR, '(A,I0)') 'Analytic test failed for task ', MYPROC
  WRITE(NERR, '(1E9.2,A3,1E9.2)') ZMAX_ERROR, ' > ', RTOLERANCE
  WRITE(NERR, '(A)') '*******************************'
  FLUSH(NERR)
  CALL ABOR1("Analytic test failed")
ELSE
  WRITE(NOUT, '(A)') '*******************************'
  WRITE(NERR, '(A,I0)') 'Analytic test passed for task ', MYPROC
  WRITE(NOUT, '(1E9.2,A3,1E9.2)') ZMAX_ERROR, ' <= ', RTOLERANCE
  WRITE(NOUT, '(A)') '*******************************'
ENDIF

!===================================================================================================

CONTAINS

!===================================================================================================

SUBROUTINE PARSE_GRID(CGRID,NDGL,NLOEN)

  CHARACTER(LEN=*) :: CGRID
  INTEGER, INTENT(INOUT) :: NDGL
  INTEGER, INTENT(INOUT), ALLOCATABLE :: NLOEN(:)
  INTEGER :: IOS
  INTEGER :: GAUSSIAN_NUMBER
  READ(CGRID(2:LEN_TRIM(CGRID)),*,IOSTAT=IOS) GAUSSIAN_NUMBER
  IF (IOS==0) THEN
    NDGL = 2 * GAUSSIAN_NUMBER
    ALLOCATE(NLOEN(NDGL))
    IF (CGRID(1:1) == 'F') THEN ! Regular Gaussian grid
      NLOEN(:) = GAUSSIAN_NUMBER * 4
      RETURN
    ENDIF
    IF (CGRID(1:1) == 'O') THEN ! Octahedral Gaussian grid
      DO I = 1, NDGL / 2
        NLOEN(I) = 20 + 4 * (I - 1)
        NLOEN(NDGL - I + 1) = NLOEN(I)
      END DO
      RETURN
    ENDIF
  ENDIF
  CALL PARSING_FAILED("ERROR: Unsupported grid specified: "// TRIM(CGRID))

END SUBROUTINE

!===================================================================================================

FUNCTION GET_REAL_VALUE(CNAME, IARG) RESULT(VALUE)

  REAL :: VALUE
  CHARACTER(LEN=*), INTENT(IN) :: CNAME
  INTEGER, INTENT(INOUT) :: IARG
  CHARACTER(LEN=128) :: CARG
  INTEGER :: STAT

  CARG = GET_STR_VALUE(CNAME, IARG)
  CALL STR2REAL(CARG, VALUE, STAT)

  IF (STAT /= 0) THEN
    CALL PARSING_FAILED("Invalid argument for " // TRIM(CNAME) // ": " // TRIM(CARG))
  END IF

END FUNCTION

!===================================================================================================

FUNCTION GET_INT_VALUE(CNAME, IARG) RESULT(VALUE)

  INTEGER :: VALUE
  CHARACTER(LEN=*), INTENT(IN) :: CNAME
  INTEGER, INTENT(INOUT) :: IARG
  CHARACTER(LEN=128) :: CARG
  INTEGER :: STAT

  CARG = GET_STR_VALUE(CNAME, IARG)
  CALL STR2INT(CARG, VALUE, STAT)

  IF (STAT /= 0) THEN
    CALL PARSING_FAILED("Invalid argument for " // TRIM(CNAME) // ": " // TRIM(CARG))
  END IF

END FUNCTION

!===================================================================================================

FUNCTION GET_STR_VALUE(CNAME, IARG) RESULT(VALUE)

  CHARACTER(LEN=128) :: VALUE
  CHARACTER(LEN=*), INTENT(IN) :: CNAME
  INTEGER, INTENT(INOUT) :: IARG

  IARG = IARG + 1
  CALL GET_COMMAND_ARGUMENT(IARG, VALUE)

  IF (VALUE == "") THEN
    CALL PARSING_FAILED("Invalid argument for " // TRIM(CNAME) // ": no value provided")
  END IF

END FUNCTION

!===================================================================================================

SUBROUTINE PARSING_FAILED(MESSAGE)

  CHARACTER(LEN=*), INTENT(IN) :: MESSAGE
  IF (LUSE_MPI) CALL MPL_INIT(LDINFO=.FALSE.)
  IF (EC_MPIRANK() == 0) THEN
    WRITE(NERR,"(A)") TRIM(MESSAGE)
    CALL PRINT_HELP(UNIT=NERR)
  ENDIF
  IF (LUSE_MPI) CALL MPL_END(LDMEMINFO=.FALSE.)
  STOP

END SUBROUTINE

!===================================================================================================

SUBROUTINE GET_COMMAND_LINE_ARGUMENTS(NSMAX, CGRID, NFLD, LUSEFLT, NPROMA, VERBOSITY, NPRTRV, &
  &                                   NPRTRW, LIMAG, RTOLERANCE)

  INTEGER, INTENT(INOUT) :: NSMAX           ! Spectral truncation
  CHARACTER(LEN=16), INTENT(INOUT) :: CGRID ! Grid
  INTEGER, INTENT(INOUT) :: NFLD            ! Number of scalar fields
  LOGICAL, INTENT(INOUT) :: LUSEFLT         ! Use fast Legendre transforms
  INTEGER, INTENT(INOUT) :: NPROMA          ! NPROMA
  INTEGER, INTENT(INOUT) :: VERBOSITY       ! Level of verbosity
  INTEGER, INTENT(INOUT) :: NPRTRV          ! Size of V set (spectral decomposition)
  INTEGER, INTENT(INOUT) :: NPRTRW          ! Size of W set (spectral decomposition)
  LOGICAL, INTENT(INOUT) :: LIMAG           ! test imaginary part
  REAL(JPRD), INTENT(INOUT) :: RTOLERANCE      ! relative error tolerance for analytic solutions

  CHARACTER(LEN=128) :: CARG          ! Storage variable for command line arguments
  INTEGER            :: IARG = 1      ! Argument index
  INTEGER            :: STAT          ! For storing success status of string->integer conversion
  INTEGER            :: MYPROC

  DO WHILE (IARG <= COMMAND_ARGUMENT_COUNT())
    CALL GET_COMMAND_ARGUMENT(IARG, CARG)

    SELECT CASE(CARG)
      ! Parse help argument
      CASE('-h', '--help')
        IF (LUSE_MPI) CALL MPL_INIT(LDINFO=.FALSE.)
        IF (EC_MPIRANK()==0) CALL PRINT_HELP()
        IF (LUSE_MPI) CALL MPL_END(LDMEMINFO=.FALSE.)
        STOP
      ! Parse verbosity argument
      CASE('-v')
        VERBOSITY = 1
      ! Parse spectral truncation argument
      CASE('-t', '--truncation')
        NSMAX = GET_INT_VALUE('-t', IARG)
        IF (NSMAX < 1) THEN
          CALL PARSING_FAILED("Invalid argument for -t: must be > 0")
        END IF
      CASE('-g', '--grid'); CGRID = GET_STR_VALUE('-g', IARG)
      CASE('-f', '--nfld'); NFLD = GET_INT_VALUE('-f', IARG)
      CASE('--imaginary'); LIMAG = .TRUE.
      CASE('--flt'); LUSEFLT = .TRUE.
      CASE('--nproma'); NPROMA = GET_INT_VALUE('--nproma', IARG)
      CASE('--nprtrv'); NPRTRV = GET_INT_VALUE('--nprtrv', IARG)
      CASE('--nprtrw'); NPRTRW = GET_INT_VALUE('--nprtrw', IARG)
      CASE('--tolerance'); RTOLERANCE = GET_REAL_VALUE('--tolerance', IARG)
      CASE DEFAULT
        CALL PARSING_FAILED("Unrecognised argument: " // TRIM(CARG))

    END SELECT
    IARG = IARG + 1
  END DO

END SUBROUTINE GET_COMMAND_LINE_ARGUMENTS

!===================================================================================================

FUNCTION CUBIC_OCTAHEDRAL_GAUSSIAN_GRID(NSMAX) RESULT(CGRID)

  CHARACTER(LEN=16) :: CGRID
  INTEGER, INTENT(IN) :: NSMAX
  WRITE(CGRID,'(A,I0)') 'O',NSMAX+1

END FUNCTION

!===================================================================================================

FUNCTION CUBIC_FULL_GRID(NSMAX) RESULT(CGRID)

  CHARACTER(LEN=16) :: CGRID
  INTEGER, INTENT(IN) :: NSMAX
  WRITE(CGRID,'(A,I0)') 'F',NSMAX+1

END FUNCTION

!===================================================================================================

SUBROUTINE STR2INT(STR, INT, STAT)

  CHARACTER(LEN=*), INTENT(IN) :: STR
  INTEGER, INTENT(OUT) :: INT
  INTEGER, INTENT(OUT) :: STAT
  READ(STR, *, IOSTAT=STAT) INT

END SUBROUTINE STR2INT

!===================================================================================================

SUBROUTINE STR2REAL(STR, REAL, STAT)

  CHARACTER(LEN=*), INTENT(IN) :: STR
  REAL, INTENT(OUT) :: REAL
  INTEGER, INTENT(OUT) :: STAT
  READ(STR, *, IOSTAT=STAT) REAL

END SUBROUTINE STR2REAL

!===================================================================================================

SUBROUTINE PRINT_HELP(UNIT)

  INTEGER, OPTIONAL :: UNIT
  INTEGER :: NOUT = 6
  IF (PRESENT(UNIT)) THEN
    NOUT = UNIT
  ENDIF

  WRITE(NOUT, "(A)") ""

  IF (JPRB == JPRD) THEN
    WRITE(NOUT, "(A)") "NAME    ectrans-benchmark-dp"
  ELSE
    WRITE(NOUT, "(A)") "NAME    ectrans-benchmark-sp"
  END IF
  WRITE(NOUT, "(A)") ""

  WRITE(NOUT, "(A)") "DESCRIPTION"
  WRITE(NOUT, "(A)") "        This program tests ecTrans by transforming fields back and forth&
    & between spectral "
  IF (JPRB == JPRD) THEN
    WRITE(NOUT, "(A)") "        space and grid-point space (double-precision version)"
  ELSE
    WRITE(NOUT, "(A)") "        space and grid-point space (single-precision version)"
  END IF
  WRITE(NOUT, "(A)") ""

  WRITE(NOUT, "(A)") "USAGE"
  IF (JPRB == JPRD) THEN
    WRITE(NOUT, "(A)") "        ectrans-benchmark-dp [options]"
  ELSE
    WRITE(NOUT, "(A)") "        ectrans-benchmark-sp [options]"
  END IF
  WRITE(NOUT, "(A)") ""

  WRITE(NOUT, "(A)") "OPTIONS"
  WRITE(NOUT, "(A)") "    -h, --help          Print this message"
  WRITE(NOUT, "(A)") "    -v                  Run with verbose output"
  WRITE(NOUT, "(A)") "    -t, --truncation T  Run with this triangular spectral truncation"
  WRITE(NOUT, "(A)") "                        (default = 79)"
  WRITE(NOUT, "(A)") "    -g, --grid GRID     Run with this grid. Possible values: O<N>, F<N>"
  WRITE(NOUT, "(A)") "                        If not specified, O<N> is used with N=truncation+1"
  WRITE(NOUT, "(A)") "                        (cubic relation)"
  WRITE(NOUT, "(A)") "    -n, --niter NITER   Run for this many inverse/direct transform"
  WRITE(NOUT, "(A)") "                        iterations (default = 10)"
  WRITE(NOUT, "(A)") "    -f, --nfld NFLD     Number of scalar fields (default = 1)"
  WRITE(NOUT, "(A)") "    --imaginary         Test imaginary part"
  WRITE(NOUT, "(A)") "    --tolerance         Test is passed if largest relative lmax-error is"
  WRITE(NOUT, "(A)") "                        smaller than this tolerance (real value)"
  WRITE(NOUT, "(A)") "    --vordiv            Also transform vorticity-divergence to wind"
  WRITE(NOUT, "(A)") "    --scders            Compute scalar derivatives (default off)"
  WRITE(NOUT, "(A)") "    --flt               Run with fast Legendre transforms (default off)"
  WRITE(NOUT, "(A)") "    --nproma NPROMA     Run with NPROMA (default no blocking: NPROMA=ngptot)"
  WRITE(NOUT, "(A)") "    --norms             Calculate and print spectral norms of transformed"
  WRITE(NOUT, "(A)") "                        fields"
  WRITE(NOUT, "(A)") "                        The computation of spectral norms will skew overall"
  WRITE(NOUT, "(A)") "                        timings"
  WRITE(NOUT, "(A)") "    --nprtrv            Size of V set in spectral decomposition"
  WRITE(NOUT, "(A)") "    --nprtrw            Size of W set in spectral decomposition"
  WRITE(NOUT, "(A)") "    -c, --check VALUE   The multiplier of the machine epsilon used as a"
  WRITE(NOUT, "(A)") "                        tolerance for correctness checking"
  WRITE(NOUT, "(A)") ""

END SUBROUTINE PRINT_HELP

!===================================================================================================

SUBROUTINE INITIALIZE_SPECTRAL_ARRAY(NSMAX, KZONAL, KTOTAL, PSPSCALAR)

  INTEGER,            INTENT(IN)  :: NSMAX          ! Spectral truncation
  INTEGER,            INTENT(IN)  :: KZONAL         ! Zonal wavenumber
  INTEGER,            INTENT(IN)  :: KTOTAL         ! Total wavenumber
  REAL(KIND=JPRB),    INTENT(OUT) :: PSPSCALAR(:,:) ! Input spectral array

  INTEGER :: JFLD
  INTEGER :: INDEX, NUM_MY_ZON_WNS
  INTEGER, ALLOCATABLE :: MY_ZON_WNS(:), NASM0(:)

  ! Get zonal wavenumbers this rank is responsible for
  CALL TRANS_INQ(KNUMP=NUM_MY_ZON_WNS)
  ALLOCATE(MY_ZON_WNS(NUM_MY_ZON_WNS))
  CALL TRANS_INQ(KMYMS=MY_ZON_WNS)

  ! First initialise all spectral coefficients to zero
  PSPSCALAR(:,:) = 0.0

  ! If rank is responsible for the chosen zonal wavenumber...
  IF (ANY(MY_ZON_WNS == KZONAL)) THEN
    ! Get array of spectral array addresses (this maps (m, n=m) to array index)
    ALLOCATE(NASM0(0:NSMAX))
    CALL TRANS_INQ(KASM0=NASM0)

    ! Find out local array index of chosen spherical harmonic
    INDEX = NASM0(KZONAL) + 2 * (KTOTAL - KZONAL)

    ! Set just that element to a constant value
    DO JFLD = 1, SIZE(PSPSCALAR, 1)
      PSPSCALAR(JFLD, INDEX) = 1.0_JPRB
    END DO
  END IF

end subroutine INITIALIZE_SPECTRAL_ARRAY

!===================================================================================================

FUNCTION DETECT_MPIRUN() RESULT(LMPI_REQUIRED)
  LOGICAL :: LMPI_REQUIRED
  INTEGER :: ILEN
  INTEGER, PARAMETER :: NVARS = 5
  CHARACTER(LEN=32), DIMENSION(NVARS) :: CMPIRUN_DETECT
  CHARACTER(LEN=4) :: CLENV_DR_HOOK_ASSERT_MPI_INITIALIZED
  INTEGER :: IVAR

  ! Environment variables that are set when mpirun, srun, aprun, ... are used
  CMPIRUN_DETECT(1) = 'OMPI_COMM_WORLD_SIZE'  ! openmpi
  CMPIRUN_DETECT(2) = 'ALPS_APP_PE'           ! cray pe
  CMPIRUN_DETECT(3) = 'PMI_SIZE'              ! intel
  CMPIRUN_DETECT(4) = 'SLURM_NTASKS'          ! slurm
  CMPIRUN_DETECT(5) = 'ECTRANS_USE_MPI'       ! forced

  LMPI_REQUIRED = .FALSE.
  DO IVAR = 1, NVARS
    CALL GET_ENVIRONMENT_VARIABLE(NAME=TRIM(CMPIRUN_DETECT(IVAR)), LENGTH=ILEN)
    IF (ILEN > 0) THEN
      LMPI_REQUIRED = .TRUE.
      EXIT ! break
    ENDIF
  ENDDO
END FUNCTION

!===================================================================================================

END PROGRAM ANALYTIC_TEST

!===================================================================================================
