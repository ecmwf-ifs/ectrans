SUBROUTINE SETUP_TRANS0(KOUT,KERR,KPRINTLEV,KMAX_RESOL,&
&                       KPRGPNS,KPRGPEW,KPRTRW,LDALLOPERM)

!**** *SETUP_TRANS* - General setup routine for transform package

!     Purpose.
!     --------
!     Resolution independent part of setup of transform package
!     Has to called BEFORE SETUP_TRANS

!**   Interface.
!     ----------
!     CALL SETUP_TRANS0(...)

!     Explicit arguments : All arguments are optional, [..] default value
!     -------------------
!     KOUT - Unit number for listing output [6]
!     KERR - Unit number for error messages [0]
!     KPRINTLEV - level of output to KOUT, 0->no output,1->normal,2->debug [0]
!     KMAX_RESOL - maximum number of different resolutions for this run [1]
!     KPRGPNS - splitting level in N-S direction in grid-point space [1]
!     KPRGPEW - splitting level in E-W direction in grid-point space [1]
!     KPRTRW  - splitting level in wave direction in spectral space [1]
!     LDALLOPERM - allocate some arrays permenately (OpenMP issue) [false]

!     The total number of (MPI)-processors has to be equal to KPRGPNS*KPRGPEW

!     Method.
!     -------

!     Externals.  SUMP_TRANS0 - initial setup routine
!     ----------  

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03

!     ------------------------------------------------------------------

#include "tsmbkind.h"

!ifndef INTERFACE

USE TPM_GEN
USE TPM_DISTR

USE SUMP_TRANS0_MOD

!endif INTERFACE

IMPLICIT NONE

INTEGER_M ,OPTIONAL,INTENT(IN) :: KOUT,KERR,KPRINTLEV,KMAX_RESOL
INTEGER_M ,OPTIONAL,INTENT(IN) :: KPRGPNS,KPRGPEW,KPRTRW
LOGICAL   ,OPTIONAL,INTENT(IN) :: LDALLOPERM

!ifndef INTERFACE

LOGICAL :: LLP1,LLP2

!     ------------------------------------------------------------------

IF(MSETUP0 /= 0) THEN
  CALL ABOR1('SETUP_TRANS0: SETUP_TRANS0 MAY ONLY BE CALLED ONCE')
ENDIF

! Default values

NOUT = 6
NERR = 0
NPRINTLEV = 0
NMAX_RESOL = 1
NPRGPNS = 1
NPRGPEW = 1
NPRTRW = 1
LALLOPERM = .FALSE.
! Optional arguments

IF(PRESENT(KOUT)) THEN
  NOUT = KOUT
ENDIF
IF(PRESENT(KERR)) THEN
  NERR = KERR
ENDIF
IF(PRESENT(KPRINTLEV)) THEN
  NPRINTLEV = KPRINTLEV
ENDIF

LLP1 = NPRINTLEV>0
LLP2 = NPRINTLEV>1
IF(LLP1) WRITE(NOUT,*) '=== ENTER ROUTINE SETUP_TRANS0 ==='

IF(PRESENT(KMAX_RESOL))THEN
  NMAX_RESOL = KMAX_RESOL
ENDIF
IF(PRESENT(KPRGPNS)) THEN
  NPRGPNS = KPRGPNS
ENDIF
IF(PRESENT(KPRGPEW)) THEN
  NPRGPEW = KPRGPEW
ENDIF
IF(PRESENT(KPRTRW)) THEN
  NPRTRW = KPRTRW
ENDIF
IF(PRESENT(LDALLOPERM)) THEN
  LALLOPERM = LDALLOPERM
ENDIF

! Initial setup
CALL SUMP_TRANS0

! Setup level 0 complete
MSETUP0 = 1

!     ------------------------------------------------------------------

!endif INTERFACE

END SUBROUTINE SETUP_TRANS0



