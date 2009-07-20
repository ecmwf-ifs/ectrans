MODULE CPLEDN_MOD
CONTAINS
SUBROUTINE CPLEDN(KN,KODD,PFN,PX,PDX,KFLAG,PW,PXN,PDXN,PXMOD)

!**** *CPLEDN* - Routine to perform a single Newton iteration step to find
!                the zero of the ordinary Legendre polynomial of degree N

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *CPLEDN(KN,KDBLE,PX,PDX,KFLAG,PW,PXN,PDXN,PXMOD)*

!        Explicit arguments :
!        --------------------
!          KN       :  Degree of the Legendre polynomial              (in)
!          KODD     :  odd or even number of latitudes                (in)
!          PFN      :  Fourier coefficients of series expansion       (in)
!                      for the ordinary Legendre polynomials
!          PX       :  abcissa where the computations are performed   (in)
!          PDX      :  id in double precision                         (in)
!          KFLAG    :  When KFLAG.EQ.1 computes the weights           (in)
!          PW       :  Weight of the quadrature at PXN                (out)
!          PXN      :  new abscissa (Newton iteration)                (out)
!          PDXN     :  id in double precision                         (out)
!          PXMOD    :  PXN-PX                                         (out)

!        Implicit arguments :
!        --------------------
!       None

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        None

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 87-10-15
!        Michel Rochas, 90-08-30 (Lobatto+cleaning)
!        K. Yessad (Sep 2008): cleaning, improve comments.
!        Nils Wedi + Mats Hamrud, 2009-02-05 revised following Swarztrauber, 2002
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KN
INTEGER(KIND=JPIM),INTENT(IN) :: KODD
REAL(KIND=JPRB),INTENT(IN)    :: PFN(0:KN/2)
REAL(KIND=JPRB),INTENT(IN)    :: PX
REAL(KIND=JPRB),INTENT(IN)    :: PDX
INTEGER(KIND=JPIM),INTENT(IN) :: KFLAG
REAL(KIND=JPRB),INTENT(OUT)   :: PW
REAL(KIND=JPRB),INTENT(OUT)   :: PXN
REAL(KIND=JPRB),INTENT(OUT)   :: PDXN
REAL(KIND=JPRB),INTENT(OUT)   :: PXMOD

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZDLX,ZDLK,ZDLKM1,ZDLKM2,ZDLLDN,ZDLXN,ZDLMOD

INTEGER(KIND=JPIM), PARAMETER :: JPKS=KIND(PX)
INTEGER(KIND=JPIM), PARAMETER :: JPKD=KIND(PDX)

INTEGER(KIND=JPIM) :: JN, IK

!      -----------------------------------------------------------------

!*       1. NEWTON ITERATION STEP.
!           ----------------------

ZDLX = PDX

ZDLK = 0.0_JPRB
IF( KODD==0 ) ZDLK=0.5_JPRB*PFN(0)
ZDLXN = 0.0_JPRB
ZDLLDN = 0.0_JPRB
IK=1

IF(KFLAG == 0)THEN
  DO JN=2-KODD,KN,2
    ! normalised ordinary Legendre polynomial == \overbar{P_n}^0
    ZDLK = ZDLK + PFN(IK)*COS(REAL(JN,JPKD)*ZDLX)
    ! normalised derivative == d/d\theta(\overbar{P_n}^0)
    ZDLLDN = ZDLLDN - PFN(IK)*REAL(JN,JPKD)*SIN(REAL(JN,JPKD)*ZDLX)
    IK=IK+1
  ENDDO
  ! Newton method
  ZDLMOD = -ZDLK/ZDLLDN
  ZDLXN = ZDLX+ZDLMOD
  PXN = REAL(ZDLXN,JPKS)
  PDXN = ZDLXN
  PXMOD = REAL(ZDLMOD,JPKS)
ENDIF

!     ------------------------------------------------------------------

!*         2.    Computes weight.
!                ----------------

IF(KFLAG == 1)THEN
  DO JN=2-KODD,KN,2
    ! normalised derivative
    ZDLLDN = ZDLLDN - PFN(IK)*REAL(JN,JPKD)*SIN(REAL(JN,JPKD)*ZDLX)
    IK=IK+1
  ENDDO
  PW = REAL(2*KN+1,JPKD)/ZDLLDN**2
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE CPLEDN
END MODULE CPLEDN_MOD
