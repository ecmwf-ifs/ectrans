SUBROUTINE SETUP_TRANS0(KOUT,KERR,KPRINTLEV,KMAX_RESOL,KPROMATR,&
&                       KPRGPNS,KPRGPEW,KPRTRW,KCOMBFLEN,LDIMP,&
&                       LDIMP_NOOLAP,LDMPOFF)

!**** *SETUP_TRANS0* - General setup routine for transform package

!     Purpose.
!     --------
!     Resolution independent part of setup of transform package
!     Has to be called BEFORE SETUP_TRANS

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
!     KCOMBFLEN - Size of communication buffer [1800000 (*8bytes) ]
!     LDIMP - use immediate message passing [false]
!     LDIMP_NOOLAP - use immediate message passing with no overlap between
!                    communications and computations [false]
!     LDMPOFF - switch off message passing [false]

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
!        R. El Khatib 03-01-24 LDMPOFF

!     ------------------------------------------------------------------

#include "tsmbkind.h"

!ifndef INTERFACE

USE TPM_GEN
USE TPM_DISTR

USE SUMP_TRANS0_MOD
USE ABORT_TRANS_MOD

!endif INTERFACE

IMPLICIT NONE

INTEGER_M ,OPTIONAL,INTENT(IN) :: KOUT,KERR,KPRINTLEV,KMAX_RESOL,KPROMATR
INTEGER_M ,OPTIONAL,INTENT(IN) :: KPRGPNS,KPRGPEW,KPRTRW,KCOMBFLEN
LOGICAL   ,OPTIONAL,INTENT(IN) :: LDIMP
LOGICAL   ,OPTIONAL,INTENT(IN) :: LDIMP_NOOLAP
LOGICAL   ,OPTIONAL,INTENT(IN) :: LDMPOFF

!ifndef INTERFACE

LOGICAL :: LLP1,LLP2

!     ------------------------------------------------------------------

IF(MSETUP0 /= 0) THEN
  CALL ABORT_TRANS('SETUP_TRANS0: SETUP_TRANS0 MAY ONLY BE CALLED ONCE')
ENDIF

! Default values

NOUT = 6
NERR = 0
NPRINTLEV = 0
NMAX_RESOL = 1
NPRGPNS = 1
NPRGPEW = 1
NPRTRW = 1
NPROMATR = 0
NCOMBFLEN = 1800000
LIMP = .FALSE.
LIMP_NOOLAP = .FALSE.
LMPOFF = .FALSE.

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
IF(PRESENT(KPROMATR))THEN
  IF(MOD(KPROMATR,2) /= 0) THEN
    CALL ABORT_TRANS('SETUP_TRANS0: KPROMATR HAS TO BE MULTIPLE OF 2')
  ENDIF
  NPROMATR = KPROMATR
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
IF(PRESENT(KCOMBFLEN)) THEN
  NCOMBFLEN = KCOMBFLEN
ENDIF
IF(PRESENT(LDIMP)) THEN
  LIMP = LDIMP
ENDIF
IF(PRESENT(LDIMP_NOOLAP)) THEN
  LIMP_NOOLAP = LDIMP_NOOLAP
ENDIF
IF(PRESENT(LDMPOFF)) THEN
  LMPOFF = LDMPOFF
ENDIF

! Initial setup
CALL SUMP_TRANS0

! Setup level 0 complete
MSETUP0 = 1

!     ------------------------------------------------------------------

!endif INTERFACE

END SUBROUTINE SETUP_TRANS0



