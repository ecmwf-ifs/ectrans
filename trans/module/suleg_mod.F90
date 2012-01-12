MODULE SULEG_MOD
CONTAINS
SUBROUTINE SULEG


USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE MPL_MODULE

USE TPM_GEN
USE TPM_DIM
USE TPM_CONSTANTS
USE TPM_DISTR
USE TPM_FIELDS
USE TPM_FLT
USE TPM_GEOMETRY

USE SUGAW_MOD
USE SUPOL_MOD
USE SUPOLF_MOD
USE TPM_POL
USE SUTRLE_MOD
USE SETUP_GEOM_MOD
USE BUTTERFLY_ALG_MOD


!**** *SULEG * - initialize the Legendre polynomials

!     Purpose.
!     --------
!           Initialize COMMON YOMLEG

!**   Interface.
!     ----------
!        *CALL* *SULEG*

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------
!              COMMON YOMLEG

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        SUGAW  (Gaussian latitudes)
!        SUPOLM (polynomials)
!        LFI routines for external IO's
!        Called by SUGEM.

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!     
!     S.L. Belousov, Tables of normalized associated Legendre Polynomials, Pergamon Press (1962)
!     P.N. Swarztrauber, On computing the points and weights for Gauss-Legendre quadrature,
!     SIAM J. Sci. Comput. Vol. 24 (3) pp. 945-954 (2002)

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 87-10-15
!        MODIFICATION : 91-04  J.M. Piriou:
!                       - Read gaussian latitudes and PNM on LFI
!                       - If file missing, computes
!                       91-04 M.Hamrud:
!                       - IO Scheme introduced
!        MODIFICATION : 91-07-03 P.Courtier suppress derivatives
!        MODIFICATION : 91-07-03 P.Courtier computes RATATH and RACTHE
!        MODIFICATION : 91-07-03 P.Courtier change upper limit (NSMAX+1)
!        MODIFICATION : 91-07-03 P.Courtier change ordering
!     Order of the PNM in the file, as in the model :
!         - increasing wave numbers m
!         - for a given m, from n=NSMAX+1 to m
!        MODIFICATION : 92-07-02 R. Bubnova: shift RATATH calculation
!                                            to SUGEM1
!        MODIFICATION : 92-12-17 P.Courtier multitask computations
!        Modified by R. EL Khatib : 93-04-02 Set-up defaults controled by LECMWF
!        MODIFICATION : 93-03-19 D.Giard : n <= NTMAX
!        K. YESSAD    : 93-05-11 : DLMU --> global array DRMU(NDGSA:NDGEN).
!                       (not stored currently on LFI files).
!        MODIFICATION : 94-02-03 R. El Khatib : subroutine SULEG2 to write out
!                       the Leg. polynomials on workfile or LFI file
!        Modification : 94-08-31 M. Tolstykh: Setup for CUD interpolation
!        Modified by K. YESSAD (MARCH 1995): Extra-latitudes computations
!                 according to value of NDGSUR and LRPOLE only.
!                 + change fancy loop numbering.
!        Modified 98-08-10 by K. YESSAD: removal of LRPOLE option.
!          - removal of LRPOLE in YOMCT0.
!          - removal of code under LRPOLE.
!        R. El Khatib 11-Apr-2007 Emulation of vectorized quadruple precision
!                                 on NEC
!        Nils Wedi + Mats Hamrud, 2009-02-05 revised following Swarztrauber, 2002
!     ------------------------------------------------------------------



IMPLICIT NONE

!     ------------------------------------------------------------------
REAL(KIND=JPRB),ALLOCATABLE :: ZPNMG(:)
REAL(KIND=JPRB),ALLOCATABLE :: ZFN(:,:)
REAL(KIND=JPRB) :: ZLRMUZ(R%NDGL)
REAL(KIND=JPRB) :: ZW(R%NDGL)

!     ------------------------------------------------------------------
REAL(KIND=JPRB) :: ZANM

REAL(KIND=JPRB) :: ZFNN

!     LOCAL 
INTEGER(KIND=JPIM) :: INM, IM , ICOUNT,&
             &JGL,  JM, JMLOC, JN, JNM, &
             &IODD, INN, INMAX, JI, IMAXN, &
             &INX, ISL, ISTART, ITHRESHOLD, INSMAX

INTEGER(KIND=JPIM) :: ILATS,ILOOP

REAL(KIND=JPRB) :: ZEPS, ZEPS_INT_DEC
REAL(KIND=JPRB),ALLOCATABLE :: ZLFPOL(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZLPOL(:)
INTEGER(KIND=JPIM) :: IMAXCOLS
INTEGER(KIND=JPIM) :: ILATSMAX,JW,JV,J
INTEGER(KIND=JPIM) :: IDGLU, ILA, ILS, IA, IS, I

LOGICAL :: LLP1,LLP2

!     ------------------------------------------------------------------

!*       0.    Some initializations.
!              ---------------------

ZEPS = 1000._JPRB*EPSILON(ZEPS)
!ZEPS_INT_DEC = EPSILON(ZEPS)
ZEPS_INT_DEC = 1.0E-7_JPRB

LLP1 = NPRINTLEV>0
LLP2 = NPRINTLEV>1
IF(LLP1) WRITE(NOUT,*) '=== ENTER ROUTINE SULEG ==='

IF( NPROC > 1 )THEN
  CALL GSTATS(796,0)
  CALL MPL_BARRIER(CDSTRING='SULEG:')
  CALL GSTATS(796,1)
ENDIF

CALL GSTATS(140,0)
CALL GSTATS(1801,0)
ALLOCATE(F%RMU(R%NDGL))
IF (LLP2) WRITE(NOUT,9) 'F%RMU     ',SIZE(F%RMU ),SHAPE(F%RMU ) 
ALLOCATE(F%RW(R%NDGL))
IF (LLP2) WRITE(NOUT,9) 'F%RW      ',SIZE(F%RW  ),SHAPE(F%RW  ) 


!*       1.0 Initialize Fourier coefficients for ordinary Legendre polynomials
!     ------------------------------------------------------------------------

ALLOCATE(ZFN(0:R%NDGL,0:R%NDGL))
IF (LLP2) WRITE(NOUT,9) 'ZFN       ',SIZE(ZFN   ),SHAPE(ZFN   )


! determines the number of stripes in butterfly NSMAX/IMAXCOLS 
! IMAXCOLS = (R%NSMAX - 1)/4 + 1
! IMAXCOLS=64 (min flops)
IMAXCOLS=128 ! (max time)
IF( R%NSMAX == 159 ) IMAXCOLS=32
IF( R%NSMAX == 63 ) IMAXCOLS=8

! the threshold of efficiency
ITHRESHOLD=900
IF( R%NSMAX <= 63 ) ITHRESHOLD=16
IF( R%NSMAX == 159 ) ITHRESHOLD=64
IF( R%NSMAX == 799 ) ITHRESHOLD=256
IF( R%NSMAX == 1279 ) ITHRESHOLD=512
IF( R%NSMAX == 2047 ) ITHRESHOLD=600
IF( R%NSMAX >= 3999 ) ITHRESHOLD=900
ITHRESHOLD=MAX(ITHRESHOLD,IMAXCOLS+1)
S%ITHRESHOLD=ITHRESHOLD

!*       3.1   Gaussian latitudes and weights
!     ---------------------------------------

CALL INI_POL(R%NTMAX+1)

IF(.NOT.D%LGRIDONLY) THEN
  ISTART=1
ELSE
  ISTART=R%NDGL
ENDIF

INMAX=R%NDGL
! Belousov, Swarztrauber use ZFN(0,0)=SQRT(2._JPRB)
! IFS normalisation chosen to be 0.5*Integral(Pnm**2) = 1
ZFN(0,0)=2._JPRB
DO JN=ISTART,R%NDGL
  ZFNN=ZFN(0,0)
  DO JGL=1,JN
    ZFNN=ZFNN*SQRT(1._JPRB-0.25_JPRB/REAL(JGL**2,JPRB))
  ENDDO

  IODD=MOD(JN,2)
  ZFN(JN,JN)=ZFNN
  DO JGL=2,JN-IODD,2
    ZFN(JN,JN-JGL)=ZFN(JN,JN-JGL+2)*REAL((JGL-1)*(2*JN-JGL+2),JPRB)/REAL(JGL*(2*JN-JGL+1),JPRB)
  ENDDO
ENDDO

! compute latitudes and weights for original Gaussian latitudes
ZANM=SQRT(REAL(2*INMAX+1,JPRB)*REAL(INMAX**2,JPRB)/REAL(2*INMAX-1,JPRB))
INN=R%NDGL
CALL SUGAW(INN,0,INMAX,ZLRMUZ(1:INN),ZW(1:INN),ZANM,ZFN)

WRITE(NOUT,*) '=== SULEG: Finished Gaussian latitudes ==='

!*       3.2   Computes related arrays

DO JGL=1,R%NDGL
  F%RW(JGL)     = ZW(JGL)
  F%RMU(JGL)    = ZLRMUZ(JGL)
ENDDO

IF(.NOT.D%LGRIDONLY) THEN

  ALLOCATE(S%FA(D%NUMP))

  ALLOCATE(F%R1MU2(R%NDGL))
  IF (LLP2) WRITE(NOUT,9) 'F%R1MU2   ',SIZE(F%R1MU2),SHAPE(F%R1MU2 ) 
  ALLOCATE(F%RACTHE(R%NDGL))
  IF (LLP2) WRITE(NOUT,9) 'F%RACTHE  ',SIZE(F%RACTHE),SHAPE(F%RACTHE ) 

  IF( S%LUSERPNM .OR. S%LKEEPRPNM ) THEN
    ALLOCATE(F%RPNM(R%NLEI3,D%NSPOLEGL))
    IF (LLP2) WRITE(NOUT,9) 'F%RPNM    ',SIZE(F%RPNM),SHAPE(F%RPNM) 
    DO JNM=1,D%NSPOLEGL
      F%RPNM(R%NLEI3,JNM) = 0.0_JPRB
    ENDDO
  ENDIF

!*       3.2   Computes related arrays

  DO JGL=1,R%NDGL
    F%R1MU2(JGL)  = 1.0_JPRB-ZLRMUZ(JGL)*ZLRMUZ(JGL)
    F%RACTHE(JGL) = 1.0_JPRB/SQRT(1.0_JPRB-ZLRMUZ(JGL)*ZLRMUZ(JGL))/REAL(RA,JPRB)
  ENDDO

!*       3.3   Working arrays

! compute the Legendre polynomials as a function of the z_k (Gaussian Latitudes)
! this may be faster than calling supolf for each m but uses extra communication
! and the parallelism is more limited ? Nils

  IF( S%LUSERPNM ) THEN

    INSMAX = R%NTMAX+1

    IF( INSMAX /= R%NDGL) THEN
      DEALLOCATE(ZFN)
      ALLOCATE(ZFN(0:INSMAX,0:INSMAX))
      ! Belousov, Swarztrauber use ZFN(0,0)=SQRT(2._JPRB)
      ! IFS normalisation chosen to be 0.5*Integral(Pnm**2) = 1
      ZFN(0,0)=2._JPRB
      DO JN=1,INSMAX
        ZFNN=ZFN(0,0)
        DO JGL=1,JN
          ZFNN=ZFNN*SQRT(1._JPRB-0.25_JPRB/REAL(JGL**2,JPRB))
        ENDDO

        IODD=MOD(JN,2)
        ZFN(JN,JN)=ZFNN
        DO JGL=2,JN-IODD,2
          ZFN(JN,JN-JGL)=ZFN(JN,JN-JGL+2)*REAL((JGL-1)*(2*JN-JGL+2),JPRB)/REAL(JGL*(2*JN-JGL+1),JPRB)
        ENDDO
      ENDDO
    ENDIF

    ALLOCATE(ZLFPOL(0:INSMAX,0:INSMAX))
    ALLOCATE(ZPNMG(R%NSPOLEG))

    ILOOP=0
    DO JGL=D%NLATLS(MYSETW,MYSETV),D%NLATLE(MYSETW,MYSETV)

      INM = 0
      CALL SUPOL(INSMAX,ZLRMUZ(JGL),ZFN,ZLFPOL)
      DO JM=0,R%NSMAX
        DO JN=INSMAX,JM,-1
          INM = INM+1
          ZPNMG(INM) = ZLFPOL(JM,JN)
        ENDDO
      ENDDO

      CALL GSTATS(1801,2)
      ILOOP = JGL-D%NLATLS(MYSETW,MYSETV)+1
      CALL SUTRLE(ZPNMG,JGL,ILOOP)
      CALL GSTATS(1801,3)

    ENDDO

    DEALLOCATE(ZLFPOL)
    IF( ALLOCATED(ZFN) ) DEALLOCATE(ZFN)

    ILATSMAX=0
    DO JW=1,NPRTRW
      DO JV=1,NPRTRV
        ILATSMAX=MAX(ILATSMAX,D%NLATLE(JW,JV)-D%NLATLS(JW,JV)+1)
      ENDDO
    ENDDO
    ILATS=D%NLATLE(MYSETW,MYSETV)-D%NLATLS(MYSETW,MYSETV)+1
    DO J=ILATS,ILATSMAX-1
      ILOOP=ILOOP+1
      CALL GSTATS(1801,2)
      CALL SUTRLE(ZPNMG,-1,ILOOP)
      CALL GSTATS(1801,3)
    ENDDO

    DEALLOCATE(ZPNMG)

    WRITE(NOUT,*) '=== SULEG: Finished RPNM ==='

  ENDIF

  CALL SETUP_GEOM

  IMAXN=R%NTMAX+1

  DO JMLOC=1,D%NUMP
    
    IM = D%MYMS(JMLOC)
    S%IMLOC = JMLOC
    
    ISL = MAX(R%NDGNH-G%NDGLU(IM)+1,1)
    IA  = 1+MOD(R%NSMAX-IM+2,2)
    IS  = 1+MOD(R%NSMAX-IM+1,2)
    
    ILA = (R%NSMAX-IM+2)/2
    ILS = (R%NSMAX-IM+3)/2
    
    IDGLU = MIN(R%NDGNH,G%NDGLU(IM))
    
    ! anti-symmetric
    IF(MOD(IMAXN-IM,2) == 0) THEN
      INMAX=IMAXN+1
    ELSE
      INMAX=IMAXN
  ENDIF

    ALLOCATE(S%FA(JMLOC)%RPNMA(IDGLU,ILA))  
    IF( .NOT. S%LUSERPNM ) THEN 
      !$OMP PARALLEL PRIVATE(JGL,ZLPOL,JI,JN)
      IF (.NOT.ALLOCATED(ZLPOL)) ALLOCATE(ZLPOL(0:INMAX))
      !$OMP DO SCHEDULE(DYNAMIC,1)
      DO JGL=1,IDGLU
        CALL SUPOLF(IM,INMAX,ZLRMUZ(ISL+JGL-1),ZLPOL(0:INMAX),KCHEAP=3)
        DO JI=1,ILA
          JN=IM+2*(JI-1)+1
          S%FA(JMLOC)%RPNMA(JGL,ILA-JI+1)=ZLPOL(JN)
        ENDDO
      ENDDO
       
      !$OMP END DO
      IF (ALLOCATED(ZLPOL)) DEALLOCATE(ZLPOL)
      !$OMP END PARALLEL

      IF( S%LKEEPRPNM) THEN
        DO J=1,ILA
          DO I=1,IDGLU
            F%RPNM(ISL+I-1,D%NPMS(IM)+IA+(J-1)*2) = S%FA(JMLOC)%RPNMA(I,J)
          ENDDO
        ENDDO
      ENDIF
    ELSE    
      DO J=1,ILA
        DO I=1,IDGLU
          S%FA(JMLOC)%RPNMA(I,J) = F%RPNM(ISL+I-1,D%NPMS(IM)+IA+(J-1)*2)
        ENDDO
      ENDDO
    ENDIF
      
    IF( ILA > ITHRESHOLD .AND. S%LUSEFLT) THEN 
      S%LSYM = .FALSE.
      INX = IDGLU
      CALL CONSTRUCT_BUTTERFLY(ZEPS_INT_DEC,IMAXCOLS,INX,ILA,S%FA(JMLOC)%RPNMA,&
       & S%FA(JMLOC)%YBUT_STRUCT_A)
      DEALLOCATE(S%FA(JMLOC)%RPNMA)
    ENDIF
    
    ! symmetric
    IF(MOD(IMAXN-IM,2) == 0) THEN
      INMAX=IMAXN
    ELSE
      INMAX=IMAXN+1
    ENDIF
    
    ALLOCATE(S%FA(JMLOC)%RPNMS(IDGLU,ILS))
    IF( .NOT. S%LUSERPNM ) THEN 

      !$OMP PARALLEL PRIVATE(JGL,ZLPOL,JI,JN)
      IF (.NOT.ALLOCATED(ZLPOL)) ALLOCATE(ZLPOL(0:INMAX))
      !$OMP DO SCHEDULE(DYNAMIC,1)
      DO JGL=1,IDGLU
        CALL SUPOLF(IM,INMAX,ZLRMUZ(ISL+JGL-1),ZLPOL(0:INMAX),KCHEAP=2)
        DO JI=1,ILS
          JN=IM+2*(JI-1)
          S%FA(JMLOC)%RPNMS(JGL,ILS-JI+1)=ZLPOL(JN)
        ENDDO
      ENDDO
      !$OMP END DO
      IF (ALLOCATED(ZLPOL)) DEALLOCATE(ZLPOL)
      !$OMP END PARALLEL
      IF( S%LKEEPRPNM) THEN
        DO J=1,ILS
          DO I=1,IDGLU
            F%RPNM(ISL+I-1,D%NPMS(IM)+IS+(J-1)*2) = S%FA(JMLOC)%RPNMS(I,J)
          ENDDO
        ENDDO
      ENDIF
    ELSE      
      DO J=1,ILS
        DO I=1,IDGLU
          S%FA(JMLOC)%RPNMS(I,J) = F%RPNM(ISL+I-1,D%NPMS(IM)+IS+(J-1)*2)
        ENDDO
      ENDDO
    ENDIF

    IF( ILS > ITHRESHOLD .AND. S%LUSEFLT) THEN
      S%LSYM = .TRUE.
      INX = IDGLU
      CALL CONSTRUCT_BUTTERFLY(ZEPS_INT_DEC,IMAXCOLS,INX,ILS,S%FA(JMLOC)%RPNMS,&
       & S%FA(JMLOC)%YBUT_STRUCT_S)
      DEALLOCATE(S%FA(JMLOC)%RPNMS)
    ENDIF

  ENDDO ! JMLOC
    
  IF( S%LUSEFLT ) THEN    
    WRITE(NOUT,*) '=== SULEG: Finished SETUP_BUTTERFLY ==='
  ENDIF


  IF( .NOT. S%LKEEPRPNM ) THEN
    IF(S%LUSERPNM) DEALLOCATE(F%RPNM)
  ENDIF

  ICOUNT = 0
  DO JMLOC=1,D%NUMP
    IM = D%MYMS(JMLOC)
    DO JN=IM,R%NTMAX+2
      ICOUNT = ICOUNT+1
    ENDDO
  ENDDO

  ALLOCATE(F%REPSNM(ICOUNT))
  IF (LLP2) WRITE(NOUT,9) 'F%REPSNM  ',SIZE(F%REPSNM ),SHAPE(F%REPSNM ) 

  ICOUNT = 0
  DO JMLOC=1,D%NUMP
    IM = D%MYMS(JMLOC)
    DO JN=IM,R%NTMAX+2
      ICOUNT = ICOUNT+1
      F%REPSNM(ICOUNT) = SQRT(REAL(JN*JN-IM*IM,JPRB)/&
       &REAL(4*JN*JN-1,JPRB))
    ENDDO
  ENDDO

  ALLOCATE(F%RN(-1:R%NTMAX+3))
  IF (LLP2) WRITE(NOUT,9) 'F%RN      ',SIZE(F%RN     ),SHAPE(F%RN     ) 
  ALLOCATE(F%RLAPIN(-1:R%NSMAX+2))
  IF (LLP2) WRITE(NOUT,9) 'F%RLAPIN  ',SIZE(F%RLAPIN ),SHAPE(F%RLAPIN ) 
  ALLOCATE(F%NLTN(-1:R%NTMAX+3))
  IF (LLP2) WRITE(NOUT,9) 'F%NLTN    ',SIZE(F%NLTN ),SHAPE(F%NLTN ) 

  DO JN=-1,R%NTMAX+3
    F%RN(JN) = REAL(JN,JPRB)
    F%NLTN(JN) = R%NTMAX+2-JN
  ENDDO
  F%RLAPIN(:)  = 0.0_JPRB
  F%RLAPIN(0)  = 0.0_JPRB
  F%RLAPIN(-1) = 0.0_JPRB
  DO JN=1,R%NSMAX+2
    F%RLAPIN(JN)=-(REAL(RA,JPRB)*REAL(RA,JPRB)/REAL(JN*(JN+1),JPRB))
  ENDDO

ENDIF
CALL GSTATS(1801,1)
CALL GSTATS(140,1)

!     ------------------------------------------------------------------
9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

CALL END_POL

END SUBROUTINE SULEG
END MODULE SULEG_MOD
