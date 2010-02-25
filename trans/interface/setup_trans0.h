SUBROUTINE SETUP_TRANS0(KOUT,KERR,KPRINTLEV,KMAX_RESOL,KPROMATR,&
&                       KPRGPNS,KPRGPEW,KPRTRW,KCOMBFLEN,&
&                       LDMPOFF,LDSYNC_TRANS,&
&                       LDEQ_REGIONS,K_REGIONS_NS,K_REGIONS_EW,K_REGIONS,&
&                       PRAD)

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
!     LDMPOFF - switch off message passing [false]
!     LDSYNC_TRANS - switch to activate barrier before transforms [false]
!     LDEQ_REGIONS - true if new eq_regions partitioning [false]
!     K_REGIONS    - Number of regions (1D or 2D partitioning)
!     K_REGIONS_NS - Maximum number of NS partitions
!     K_REGIONS_EW - Maximum number of EW partitions
!     PRAD         - Radius of the planet

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
!        G. Mozdzynski 2006-09-13 LDEQ_REGIONS
!        N. Wedi  2009-11-30 add radius

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN) :: KOUT,KERR,KPRINTLEV,KMAX_RESOL,KPROMATR
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN) :: KPRGPNS,KPRGPEW,KPRTRW,KCOMBFLEN
LOGICAL   ,OPTIONAL,INTENT(IN) :: LDMPOFF
LOGICAL   ,OPTIONAL,INTENT(IN) :: LDSYNC_TRANS
LOGICAL   ,OPTIONAL,INTENT(IN) :: LDEQ_REGIONS
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(OUT) :: K_REGIONS(:)
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(OUT) :: K_REGIONS_NS
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(OUT) :: K_REGIONS_EW
REAL(KIND=JPRB),OPTIONAL,INTENT(IN) :: PRAD

END SUBROUTINE SETUP_TRANS0



