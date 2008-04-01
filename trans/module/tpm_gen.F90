MODULE TPM_GEN

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

INTEGER(KIND=JPIM) :: NOUT            ! Unit number for "standard" output
INTEGER(KIND=JPIM) :: NERR            ! Unit number for error messages
INTEGER(KIND=JPIM) :: NPRINTLEV       ! Printing level, 0=no print, 1=standard,2=debug

INTEGER(KIND=JPIM) :: MSETUP0 = 0     ! Control of setup calls
INTEGER(KIND=JPIM) :: NMAX_RESOL = 0  ! Maximum allowed number of resolutions
INTEGER(KIND=JPIM) :: NCUR_RESOL = 0  ! Current resolution
INTEGER(KIND=JPIM) :: NDEF_RESOL = 0  ! Number of defined resolutions
INTEGER(KIND=JPIM) :: NPROMATR        ! Packet size for transform (in no of fields)
                             ! NPROMATR=0 means do all fields together (dflt)

LOGICAL   :: LALLOPERM       ! Allocate some shared data structures permanently
LOGICAL   :: LMPOFF          ! true: switch off message passing
LOGICAL   :: LSYNC_TRANS     ! true: activate barriers in trmtol and trltom

END MODULE TPM_GEN
