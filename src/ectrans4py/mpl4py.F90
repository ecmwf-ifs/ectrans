SUBROUTINE MPL_INIT4PY(KRANK, KSIZE)
! ** PURPOSE
!    Initialise FIAT's MPL (Message Passing Library) for distributed ecTrans.
!    MPL_INIT detects an already-initialised MPI (e.g. from mpi4py) and attaches
!    to MPI_COMM_WORLD, so this must be called *after* MPI is up.
!
! ** RETURNS
!    KRANK : this task's MPL rank (1-based)
!    KSIZE : total number of MPI tasks
USE ISO_FORTRAN_ENV, ONLY: INT64
USE MPL_MODULE, ONLY: MPL_INIT, MPL_MYRANK, MPL_NPROC
IMPLICIT NONE
INTEGER(KIND=INT64), INTENT(OUT) :: KRANK, KSIZE
CALL MPL_INIT(LDINFO=.FALSE.)
KRANK = MPL_MYRANK()
KSIZE = MPL_NPROC()
END SUBROUTINE MPL_INIT4PY


SUBROUTINE MPL_END4PY()
! ** PURPOSE
!    Finalise FIAT's MPL. (Does not finalise MPI itself.)
USE MPL_MODULE, ONLY: MPL_END
IMPLICIT NONE
CALL MPL_END(LDMEMINFO=.FALSE.)
END SUBROUTINE MPL_END4PY


SUBROUTINE MPL_INIT_COMM4PY(KCOMM, KRANK, KSIZE)
! ** PURPOSE
!    Initialise FIAT's MPL like MPL_INIT4PY, but on the user communicator KCOMM instead of
!    MPI_COMM_WORLD (FIAT's LMPLUSERCOMM and MPLUSERCOMM, as transi's trans_set_mpi_comm
!    uses): every MPL communication, and so all of ecTrans, then stays within the tasks of
!    KCOMM, and the other tasks of MPI_COMM_WORLD take no part. Call it instead of
!    MPL_INIT4PY, after MPI is up, on the tasks of KCOMM. MPL is initialised once per
!    process: if it already is, it is left as it is.
!
! ** ARGUMENTS
!    KCOMM : Fortran handle of the communicator (mpi4py: comm.py2f())
!
! ** RETURNS
!    KRANK : this task's MPL rank (1-based)
!    KSIZE : number of MPL tasks, those of KCOMM
USE ISO_FORTRAN_ENV, ONLY: INT64
USE PARKIND1, ONLY: JPIM
USE MPL_MODULE, ONLY: MPL_INIT, MPL_MYRANK, MPL_NPROC, LMPLUSERCOMM, MPLUSERCOMM
IMPLICIT NONE
INTEGER(KIND=INT64), INTENT(IN)  :: KCOMM
INTEGER(KIND=INT64), INTENT(OUT) :: KRANK, KSIZE
LMPLUSERCOMM = .TRUE.
MPLUSERCOMM = INT(KCOMM, JPIM)
CALL MPL_INIT(LDINFO=.FALSE.)
KRANK = MPL_MYRANK()
KSIZE = MPL_NPROC()
END SUBROUTINE MPL_INIT_COMM4PY
