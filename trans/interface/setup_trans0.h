SUBROUTINE SETUP_TRANS0(KOUT,KERR,KPRINTLEV,KMAX_RESOL,KPROMATR,&
&                       KPRGPNS,KPRGPEW,KPRTRW,KCOMBFLEN,LDALLOPERM,LDIMP,&
&                       LDIMP_NOOLAP)

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
!     LDALLOPERM - allocate some arrays permenately (OpenMP issue) [false]
!     LDIMP - use immediate message passing [false]
!     LDIMP_NOOLAP - use immediate message passing with no overlap between
!                    communications and computations [false]

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


IMPLICIT NONE

INTEGER_M ,OPTIONAL,INTENT(IN) :: KOUT,KERR,KPRINTLEV,KMAX_RESOL,KPROMATR
INTEGER_M ,OPTIONAL,INTENT(IN) :: KPRGPNS,KPRGPEW,KPRTRW,KCOMBFLEN
LOGICAL   ,OPTIONAL,INTENT(IN) :: LDALLOPERM
LOGICAL   ,OPTIONAL,INTENT(IN) :: LDIMP
LOGICAL   ,OPTIONAL,INTENT(IN) :: LDIMP_NOOLAP


END SUBROUTINE SETUP_TRANS0



