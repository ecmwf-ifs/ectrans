MODULE TPM_GEN

#include "tsmbkind.h"

IMPLICIT NONE

SAVE

INTEGER_M :: NOUT            ! Unit number for "standard" output
INTEGER_M :: NERR            ! Unit number for error messages
INTEGER_M :: NPRINTLEV       ! Printing level, 0=no print, 1=standard,2=debug

INTEGER_M :: MSETUP0 = 0     ! Control of setup calls
INTEGER_M :: NMAX_RESOL = 0  ! Maximum allowed number of resolutions
INTEGER_M :: NCUR_RESOL = 0  ! Current resolution
INTEGER_M :: NDEF_RESOL = 0  ! Number of defined resolutions
INTEGER_M :: NPROMATR        ! Packet size for transform (in no of fields)
                             ! NPROMATR=0 means do all fields together (dflt)

LOGICAL   :: LALLOPERM       ! Allocate some shared data structures permanently
LOGICAL   :: LIMP            ! true: use immediate message passing 

END MODULE TPM_GEN
