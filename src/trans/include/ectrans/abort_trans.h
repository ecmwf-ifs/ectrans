INTERFACE
MODULE SUBROUTINE ABORT_TRANS(CDTEXT) 
USE TPM_GEN  , ONLY : NOUT,NERR
USE TPM_DISTR, ONLY : NPROC,MYPROC
USE MPL_MODULE, ONLY : MPL_ABORT
USE SDL_MOD, ONLY : SDL_TRACEBACK, SDL_SRLABORT

IMPLICIT NONE

CHARACTER(LEN=*),INTENT(IN) :: CDTEXT
END SUBROUTINE ABORT_TRANS
END INTERFACE 
