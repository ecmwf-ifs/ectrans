SUBROUTINE SPECNORM(PSPEC,KVSET,KMASTER,KRESOL,PMET,PNORM)

!**** *SPECNORM* - Compute global spectral norms

!     Purpose.
!     --------
!        Interface routine for computing spectral norms

!**   Interface.
!     ----------
!     CALL SPECNORM(...)

!     Explicit arguments : All arguments optional
!     -------------------- 
!     PSPEC(:,:)  - Spectral array
!     KVSET(:)    - "B-Set" for each field
!     KMASTER     - processor to recieve norms
!     KRESOL      - resolution tag  which is required ,default is the
!                   first defined resulution (input)
!     PMET(:)     - metric
!     PNORM(:)    - Norms (output for processor KMASTER)
!
!     ------------------------------------------------------------------

#include "tsmbkind.h"

!ifndef INTERFACE

USE TPM_GEN
USE TPM_DIM
USE TPM_DISTR

USE SET_RESOL_MOD
USE SPNORM_CTL_MOD

!endif INTERFACE

IMPLICIT NONE

! Declaration of arguments


REAL_B    ,OPTIONAL, INTENT(IN)  :: PSPEC(:,:)
INTEGER_M ,OPTIONAL, INTENT(IN)  :: KVSET(:)
INTEGER_M ,OPTIONAL, INTENT(IN)  :: KMASTER
REAL_B    ,OPTIONAL, INTENT(IN)  :: PMET(:)
REAL_B    ,OPTIONAL, INTENT(OUT) :: PNORM(:)
INTEGER_M ,OPTIONAL, INTENT(IN)  :: KRESOL
!ifndef INTERFACE

INTEGER_M :: IMASTER,IFLD,IFLD_G,J

!     ------------------------------------------------------------------

! Set current resolution
CALL SET_RESOL(KRESOL)

! Set defaults
IMASTER = 1
IFLD    = 0


IF(PRESENT(KMASTER)) THEN
  IMASTER = KMASTER
ENDIF

IF(PRESENT(KVSET)) THEN
  IFLD_G = UBOUND(KVSET,1)
  DO J=1,IFLD_G
    IF(KVSET(J) > NPRTRV) THEN
      WRITE(NERR,*) 'SPECNORM:KVSET(J) > NPRTRV ',J,KVSET(J),NPRTRV
      CALL ABOR1('SPECNORM:KVSET TOO LONG OR CONTAINS VALUES OUTSIDE RANGE')
    ENDIF
    IF(KVSET(J) == MYSETV) THEN
      IFLD = IFLD+1
    ENDIF
  ENDDO
ELSE
  IF(PRESENT(PSPEC)) THEN
    IFLD = UBOUND(PSPEC,1)
  ENDIF
  IFLD_G = IFLD
ENDIF

IF(NPRTRV >1) THEN
  IF(IFLD > 0 .AND. .NOT. PRESENT(KVSET)) THEN
    WRITE(NERR,*)'NPRTRV >1 AND IFLD > 0 AND NOT PRESENT(KVSET)',&
                 &NPRTRV,IFLD
    CALL ABOR1('SPECNORM: SPECIFY VERTICAL SPECTRAL DISTRIBUTION!')
  ENDIF
ENDIF
IF(MYPROC == IMASTER) THEN
  IF(.NOT. PRESENT(PNORM)) THEN
    CALL ABOR1('SPECNORM: PNORM NOT PRESENT')
  ENDIF
  IF(UBOUND(PNORM,1) < IFLD_G) THEN
    CALL ABOR1('SPECNORM: PNORM TOO SMALL')
  ENDIF
ENDIF
IF(IFLD > 0 ) THEN
  IF(.NOT. PRESENT(PSPEC)) THEN
    CALL ABOR1('SPECNORM: PSPEC NOT PRESENT')
  ENDIF
  IF(UBOUND(PSPEC,1) < IFLD) THEN
    CALL ABOR1('SPECNORM: FIRST DIMENSION OF PSPEC TOO SMALL')
  ENDIF
  IF(UBOUND(PSPEC,2) < D%NSPEC2) THEN
    CALL ABOR1('SPECNORM: FIRST DIMENSION OF PSPEC TOO SMALL')
  ENDIF
ENDIF  

CALL SPNORM_CTL(PSPEC,IFLD,IFLD_G,KVSET,IMASTER,PMET,PNORM)

!endif INTERFACE

!     ------------------------------------------------------------------

END SUBROUTINE SPECNORM

