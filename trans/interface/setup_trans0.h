SUBROUTINE SETUP_TRANS0(KOUT,KERR,KPRINTLEV,KMAX_RESOL,&
&                       KPRGPNS,KPRGPEW,KPRTRW)

!**** *SETUP_TRANS* - General setup routine for transform package

!     Purpose.
!     --------
!     Resolution independent part of setup of transform package
!     Has to called BEFORE SETUP_TRANS

!**   Interface.
!     ----------
!     CALL SETUP_TRANS0(...)

!     Explicit arguments : All arguments are optional
!     -------------------
!     KOUT - Unit number for listing output
!     KERR - Unit number for error messages
!     KPRINTLEV - level of output to listing, 0->no output,1->normal,2->debug 
!     KMAX_RESOL - maximum number of different resolutions for this run
!     KPRGPNS - splitting level in N-S direction in grid-point space
!     KPRGPEW - splitting level in E-W direction in grid-point space
!     KPRTRW  - splitting level in wave direction in spectral space

!     The total number of (MPI)-processors has to be equal to KPRGPNS*KPRGPEW

!     ------------------------------------------------------------------

#include "tsmbkind.h"


IMPLICIT NONE

INTEGER_M ,OPTIONAL,INTENT(IN) :: KOUT,KERR,KPRINTLEV,KMAX_RESOL
INTEGER_M ,OPTIONAL,INTENT(IN) :: KPRGPNS,KPRGPEW,KPRTRW


END SUBROUTINE SETUP_TRANS0



