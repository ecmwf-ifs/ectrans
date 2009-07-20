MODULE SUGAW_MOD
CONTAINS
SUBROUTINE SUGAW(KN,PFN,PL,PDL,PW)

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_CONSTANTS ,ONLY : RA

USE TPM_GEN
USE GAWL_MOD
USE ABORT_TRANS_MOD

!**** *SUGAW * - Routine to initialize the Gaussian
!                 abcissa and the associated weights

!     Purpose.
!     --------
!           Initialize arrays PL,PDL and PW (quadrature abscissas and weights)
!**   Interface.
!     ----------
!        *CALL* *SUGAW(KN,PFN,PL,PDL,PW) *

!        Explicit arguments :
!        --------------------
!           INPUT:
!              KN       :  Number of Gauss  abscissas
!             PFN       :  Fourier coefficients of series expansion for
!                          the ordinary Legendre polynomials
!           OUTPUT:
!              PL (KN)  :  abscissas of Gauss
!              PDL(KN)  :  idem in double precision
!              PW (KN)  :  Weights of the Gaussian integration

!     PL (i) is the abscissa i starting from the northern pole, it is
! the cosine of the colatitude of the corresponding row of the collocation
! grid.

!        Implicit arguments :
!        --------------------
!       None

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------

!     S.L. Belousov, Tables of normalized associated Legendre Polynomials, Pergamon Press (1962)
!     P.N. Swarztrauber, On computing the points and weights for Gauss-Legendre quadrature,
!     SIAM J. Sci. Comput. Vol. 24 (3) pp. 945-954 (2002)

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original          : 87-10-15
!        Michel Rochas     : 90-08-30
!        Philippe Courtier : 92-12-19 Multitasking
!        Ryad El Khatib    : 94-04-20 Remove unused comdecks pardim and yomdim
!        Mats Hamrud       : 94-08-12 Printing level
!        K. Yessad (Sep 2008): cleaning, improve comments.
!        Nils Wedi + Mats Hamrud, 2009-02-05 revised following Swarztrauber, 2002
!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KN
REAL(KIND=JPRB),INTENT(IN) :: PFN(0:KN/2)
REAL(KIND=JPRB),INTENT(OUT) :: PL(:)
REAL(KIND=JPRB),INTENT(OUT) :: PDL(:)
REAL(KIND=JPRB),INTENT(OUT) :: PW(:)

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZLI(KN),ZT(KN)
REAL(KIND=JPRB) :: ZREG(KN),ZMOD(KN),ZM(KN),ZR(KN)
INTEGER(KIND=JPIM) :: ITER(KN)
INTEGER(KIND=JPIM) :: IALLOW, INS2, ISYM, JGL
REAL(KIND=JPRB) :: Z, ZEPS, ZPI
LOGICAL :: LLP2

!     ------------------------------------------------------------------

!*       1. Initialization.
!           ---------------

LLP2 = .FALSE.
ZPI  = 2.0_JPRB*ASIN(1.0_JPRB)
INS2 = KN/2+MOD(KN,2)

!*       1.1 Find first approximation of the roots of the
!           Legendre polynomial of degree KN.
DO JGL=1,INS2
  Z = REAL(4*JGL-1,JPRB)*ZPI/REAL(4*KN+2,JPRB)
  PL(JGL) = Z+1.0_JPRB/(TAN(Z)*REAL(8*KN**2,JPRB))
  ZREG(JGL) = COS(Z)
  ZLI(JGL) = COS(PL(JGL))
ENDDO

!     ------------------------------------------------------------------

!*      2. Computes roots and weights for transformed theta 
!          ------------------------------------------------

ZEPS = EPSILON(Z)
CALL GSTATS(1650,0)
!!!!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JGL)
DO JGL=INS2,1,-1
  CALL GAWL(PFN,PL(JGL),PDL(JGL),PW(JGL),ZEPS,KN,ITER(JGL),ZMOD(JGL))
ENDDO
!!!!!$OMP END PARALLEL DO
CALL GSTATS(1650,1)

! convert to physical latitude space PMU
!DIR$ IVDEP
!OCL NOVREC
DO JGL=1,INS2
  PL(JGL) = COS(PL(JGL))
  PDL(JGL) = COS(PDL(JGL))
ENDDO

!DIR$ IVDEP
!OCL NOVREC
DO JGL=1,KN/2
  ISYM = KN-JGL+1
  PL(ISYM) = -PL(JGL)
  PDL(ISYM) = -PDL(JGL)
  PW(ISYM) = PW(JGL)
ENDDO

!     ------------------------------------------------------------------

!*      3. Diagnostics.
!          ------------

IF(LLP2)THEN
  DO JGL=1,INS2
    ZM(JGL) = (ACOS(PL(JGL))-ACOS(ZLI(JGL)))*RA
    ZR(JGL) = (ACOS(PL(JGL))-ACOS(ZREG(JGL)))*RA
    ZT(JGL) = ACOS(PL(JGL))*180._JPRB/ZPI
  ENDDO
ENDIF

IALLOW = 10
DO JGL=1,INS2
  IF(ITER(JGL) > IALLOW)THEN
    WRITE(UNIT=NOUT,FMT='('' CONVERGENCE FAILED IN SUGAW '')')
    WRITE(UNIT=NOUT,FMT='('' ALLOWED : '',I4,''&
     &NECESSARY : '',&
     &I4)')IALLOW,ITER(JGL)
    CALL ABORT_TRANS(' FAILURE IN SUGAW ')
  ENDIF

  IF(LLP2)THEN
    WRITE(UNIT=NOUT,FMT=&
     &'('' ROW ='',I4,'' ITERATIONS='',I4,'' ROOT='',F30.20,&
     &'' WEIGHT='',F30.20,'' MODIF :'',E8.2)')JGL,ITER(JGL),PL(JGL)&
     &,PW(JGL),PL(JGL)-ZLI(JGL)
    WRITE(UNIT=NOUT,FMT=&
     &'(10X,'' LAST INC. : '',E8.2,'' MODIF IN M : '',F10.3,&
     &'' FROM THE REGULAR GRID : '',F10.3,'' COLAT '',F10.3)')&
     &ZMOD(JGL),ZM(JGL),ZR(JGL),ZT(JGL)
  ENDIF
ENDDO

!     ------------------------------------------------------------------

END SUBROUTINE SUGAW
END MODULE SUGAW_MOD
