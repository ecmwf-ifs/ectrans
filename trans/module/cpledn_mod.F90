MODULE CPLEDN_MOD
CONTAINS
SUBROUTINE CPLEDN(KN,KODD,PFN,PX,DDX,KFLAG,PW,PXN,DDXN,PXMOD)

!**** *CPLEDN* - Routine to perform a single Newton iteration step to find
!                the zero of the ordinary Legendre polynomial of degree N

!     Purpose.
!     --------
!           

!**   Interface.
!     ----------
!        *CALL* *CPLEDN(KN,KODD,PFN,PX,DDX,KFLAG,PW,PXN,DDXN,PXMOD)*

!        Explicit arguments :
!        --------------------
!              KN       :  Degree of the Legendre polynomial
!              KODD     :  odd or even number of latitudes          
!              PFN      :  Fourier coefficients of series expansion 
!                          for the ordinary Legendre polynomials
!              PX       :  abcissa where the computations are performed
!              DDX      :  id in double precision
!              KFLAG    :  When KFLAG.EQ.1 computes the weights
!              PW       :  Weight of the quadrature at PXN
!              PXN      :  new abscissa (Newton iteration)
!              DDXN     :  id in double precision
!              PXMOD    :  PXN-PX

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
!        Nils Wedi + Mats Hamrud, 2009-02-05 revised following Swarztrauber, 2002
!     ------------------------------------------------------------------



USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER(KIND=JPIM) :: KFLAG
INTEGER(KIND=JPIM) :: KN
INTEGER(KIND=JPIM) :: KODD

!     DUMMY REAL SCALARS
REAL(KIND=JPRB) :: PW
REAL(KIND=JPRB) :: PX
REAL(KIND=JPRB) :: PXMOD
REAL(KIND=JPRB) :: PXN
REAL(KIND=JPRB) :: PFN(0:KN/2)

REAL(KIND=JPRB) :: DDX,DDXN,DLX,DLK,DLKM1,DLKM2,DLLDN,DLXN,DLMOD
REAL(KIND=JPRB) :: DLG,DLGDN

INTEGER(KIND=JPIM), PARAMETER :: JPKS=KIND(PX)
INTEGER(KIND=JPIM), PARAMETER :: JPKD=KIND(DDX)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: JN, IK

!      -----------------------------------------------------------------

DLX = DDX
DLK = 0.0_JPRB
IF( KODD==0 ) DLK=0.5_JPRB*PFN(0)
DLXN = 0.0_JPRB
DLLDN = 0.0_JPRB

!*       2.1   NEWTON ITERATION STEP.

IK=1
IF(KFLAG == 0)THEN
  DO JN=2-KODD,KN,2
! normalised ordinary Legendre polynomial == \overbar{P_n}^0
    DLK = DLK + PFN(IK)*COS(REAL(JN,JPKD)*DLX)
! normalised derivative == d/d\theta(\overbar{P_n}^0)
    DLLDN = DLLDN - PFN(IK)*REAL(JN,JPKD)*SIN(REAL(JN,JPKD)*DLX)
    IK=IK+1
  ENDDO
! Newton method
  DLMOD = -DLK/DLLDN
  DLXN = DLX+DLMOD
  PXN = REAL(DLXN,JPKS)
  DDXN = DLXN
  PXMOD = REAL(DLMOD,JPKS)
ELSE
!     ------------------------------------------------------------------

!*         3.    Computes weight.
!                ----------------
  DO JN=2-KODD,KN,2
! normalised derivative
    DLLDN = DLLDN - PFN(IK)*REAL(JN,JPKD)*SIN(REAL(JN,JPKD)*DLX)
    IK=IK+1
  ENDDO
  PW = REAL(2*KN+1,JPKD)/DLLDN**2
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE CPLEDN
END MODULE CPLEDN_MOD
