SUBROUTINE MPL_INIT4PY(KCOMM, KRANK, KSIZE)
! ** PURPOSE
!    Initialise FIAT's MPL on the user communicator KCOMM, unless KCOMM = -1, in which case
!    MPI_COMM_WORLD is used. Every MPL communication, and so all of ecTrans, stays within the tasks
!    of KCOMM, and the other tasks of MPI_COMM_WORLD take no part. MPL is initialised once per
!    process: if it already is, it is left as it is.
!
! ** ARGUMENTS
!    KCOMM : Fortran handle of the communicator (mpi4py: comm.py2f()). Pass -1 to get MPI_COMM_WORLD
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
IF (KCOMM /= -1) THEN
   LMPLUSERCOMM = .TRUE.
   MPLUSERCOMM = INT(KCOMM, JPIM)
END IF
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
