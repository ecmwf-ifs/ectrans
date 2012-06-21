MODULE SULEG_MOD
CONTAINS
SUBROUTINE SULEG

USE PARKIND1  ,ONLY : JPIM     ,JPRB     ,JPIB
USE PARKIND2  ,ONLY : JPRH

USE TPM_GEN
USE TPM_DIM
USE TPM_CONSTANTS
USE TPM_DISTR
USE TPM_FIELDS

USE SUGAW_MOD
USE SUPOL_MOD
USE SUTRLE_MOD

#if defined(NECSX) && defined(REALHUGE)
USE quad_emu
#endif

#ifdef DOC

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
!     ------------------------------------------------------------------

#endif


IMPLICIT NONE


!     ------------------------------------------------------------------
REAL(KIND=JPRB),ALLOCATABLE :: ZPNMG(:,:)

REAL(KIND=JPRH) :: DLRMU(R%NDGL)
REAL(KIND=JPRH) :: DLC(0:R%NTMAX+1,0:R%NTMAX+1)
REAL(KIND=JPRH) :: DLD(0:R%NTMAX+1,0:R%NTMAX+1)
REAL(KIND=JPRH) :: DLE(0:R%NTMAX+1,0:R%NTMAX+1)
REAL(KIND=JPRH) :: DLA(0:R%NTMAX+1),DLB(0:R%NTMAX+1),DLF(0:R%NTMAX+1)
REAL(KIND=JPRH) :: DLG(0:R%NTMAX+1),DLH(0:R%NTMAX+1),DLI(0:R%NTMAX+1)
REAL(KIND=JPRH) :: DLPOL(0:R%NTMAX+1,0:R%NTMAX+1)
!     ------------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: JPKS=KIND(ZPNMG)
INTEGER(KIND=JPIM), PARAMETER :: JPKD=KIND(DLG)

!     ------------------------------------------------------------------
REAL(KIND=JPRH) :: DA,DC,DD,DE
INTEGER(KIND=JPIM) :: KKN, KKM

!     LOCAL 
INTEGER(KIND=JPIM) :: IGLLOC, INM, IM , ICOUNT,&
             &JGL,  JM, JMLOC, JN, JNM


LOGICAL :: LLP1,LLP2

#if defined(NECSX) && defined(REALHUGE)
INTEGER(KIND=JPIB) :: IDCN(0:R%NTMAX+1,0:R%NTMAX+1)
INTEGER(KIND=JPIB) :: IDCD(0:R%NTMAX+1,0:R%NTMAX+1)
INTEGER(KIND=JPIB) :: IDDN(0:R%NTMAX+1,0:R%NTMAX+1)
INTEGER(KIND=JPIB) :: IDDD(0:R%NTMAX+1,0:R%NTMAX+1)
INTEGER(KIND=JPIB) :: IDEN(0:R%NTMAX+1,0:R%NTMAX+1)
INTEGER(KIND=JPIB) :: IDED(0:R%NTMAX+1,0:R%NTMAX+1)
#endif

DC(KKN,KKM)=SQRT( (REAL(2*KKN+1,JPKD)*REAL(KKN+KKM-1,JPKD)&
                   &*REAL(KKN+KKM-3,JPKD))&
                &/ (REAL(2*KKN-3,JPKD)*REAL(KKN+KKM,JPKD)&
                   &*REAL(KKN+KKM-2,JPKD)) )
DD(KKN,KKM)=SQRT( (REAL(2*KKN+1,JPKD)*REAL(KKN+KKM-1,JPKD)&
                   &*REAL(KKN-KKM+1,JPKD))&
                &/ (REAL(2*KKN-1,JPKD)*REAL(KKN+KKM,JPKD)&
                  &*REAL(KKN+KKM-2,JPKD)) )
DE(KKN,KKM)=SQRT( (REAL(2*KKN+1,JPKD)*REAL(KKN-KKM,JPKD))&
                &/ (REAL(2*KKN-1,JPKD)*REAL(KKN+KKM,JPKD)) )
DA(KKN,KKM)=SQRT( (REAL(2*KKN+1,JPKD)*REAL(KKN-KKM,JPKD)&
                   &*REAL(KKN+KKM,JPKD))&
                &/  REAL(2*KKN-1,JPKD) )

!     ------------------------------------------------------------------
ALLOCATE(ZPNMG(R%NSPOLEG,D%NLEI3D))

!*       0.    Some initializations.
!              ---------------------

LLP1 = NPRINTLEV>0
LLP2 = NPRINTLEV>1
IF(LLP1) WRITE(NOUT,*) '=== ENTER ROUTINE SULEG ==='

CALL GSTATS(140,0)
ALLOCATE(F%RPNM(R%NLEI3,D%NSPOLEGL))
IF (LLP2) WRITE(NOUT,9) 'F%RPNM    ',SIZE(F%RPNM),SHAPE(F%RPNM) 
ALLOCATE(F%RMU(R%NDGL))
IF (LLP2) WRITE(NOUT,9) 'F%RMU     ',SIZE(F%RMU ),SHAPE(F%RMU ) 
ALLOCATE(F%RW(R%NDGL))
IF (LLP2) WRITE(NOUT,9) 'F%RW      ',SIZE(F%RW  ),SHAPE(F%RW  ) 
ALLOCATE(F%R1MU2(R%NDGL))
IF (LLP2) WRITE(NOUT,9) 'F%R1MU2   ',SIZE(F%R1MU2),SHAPE(F%R1MU2 ) 
ALLOCATE(F%RACTHE(R%NDGL))
IF (LLP2) WRITE(NOUT,9) 'F%RACTHE  ',SIZE(F%RACTHE),SHAPE(F%RACTHE ) 

CALL GSTATS(1801,0)
DO JNM=1,D%NSPOLEGL
  F%RPNM(R%NLEI3,JNM) = 0.0_JPRB
ENDDO
CALL GSTATS(1801,1)

!     ------------------------------------------------------------------

!*       3.1   Gaussian latitudes and weights
CALL SUGAW(R%NDGL,F%RMU,DLRMU,F%RW)

!*       3.2   Computes related arrays

DO JGL=1,R%NDGL
  F%R1MU2(JGL) = REAL(1.0_JPRB-DLRMU(JGL)*DLRMU(JGL),JPKS)
  F%RACTHE(JGL) = REAL(1.0_JPRB/SQRT(1.0_JPRB-DLRMU(JGL)*DLRMU(JGL))/&
   &REAL(RA,JPKD),JPKS)
ENDDO

!*       3.3   Working arrays

#if defined(NECSX) && defined(REALHUGE)

DO JN=3,R%NTMAX+1
  DO JM=2,JN-1
    IDCN(JM,JN)=(2*JN+1)*(JN+JM-1)*(JN+JM-3)
    IDCD(JM,JN)=(2*JN-3)*(JN+JM  )*(JN+JM-2)
    IDDN(JM,JN)=(2*JN+1)*(JN+JM-1)*(JN-JM+1)
    IDDD(JM,JN)=(2*JN-1)*(JN+JM  )*(JN+JM-2)
    IDEN(JM,JN)=(2*JN+1)*(JN-JM)
    IDED(JM,JN)=(2*JN-1)*(JN+JM)
  ENDDO
ENDDO

DO JN=3,R%NTMAX+1
  CALL SQRTIVDV (IDCN(2,JN), IDCD(2,JN), JN-2, DLC(2,JN))
  CALL SQRTIVDV (IDDN(2,JN), IDDD(2,JN), JN-2, DLD(2,JN))
  CALL SQRTIVDV (IDEN(2,JN), IDED(2,JN), JN-2, DLE(2,JN))
ENDDO

#else

DO JN=3,R%NTMAX+1
  DO JM=2,JN-1
    DLC(JM,JN) = DC(JN,JM)
    DLD(JM,JN) = DD(JN,JM)
    DLE(JM,JN) = DE(JN,JM)
  ENDDO
ENDDO

#endif

DO JN=1,R%NTMAX+1
  DLA(JN) = SQRT(REAL(2*JN+1,JPKD))
  DLB(JN) = SQRT(REAL(2*JN+1,JPKD)/REAL(JN*(JN+1),JPKD))
  DLF(JN) = REAL(2*JN-1,JPKD)/REAL(JN,JPKD)
  DLG(JN) = REAL(JN-1,JPKD)/REAL(JN,JPKD)
  DLH(JN) = SQRT(REAL(2*JN+1,JPKD)/REAL(2*JN,JPKD))
  DLI(JN) = REAL(JN,JPKD)
ENDDO

CALL GSTATS(1801,0)
DO JGL=D%NLATLS(MYSETW),D%NLATLE(MYSETW)
  DLPOL(:,:) = 0.0_JPRB
  CALL SUPOL(R%NTMAX+1,DLRMU(JGL),DLPOL,DLA,DLB,DLC,DLD,DLE,DLF,DLG,DLH,DLI)
  INM = 0
  IGLLOC = JGL - D%NLATLS(MYSETW) + 1
  DO JM=0,R%NSMAX
    DO JN=R%NTMAX+1,JM,-1
      INM = INM+1
      ZPNMG(INM,IGLLOC) = REAL(DLPOL(JM,JN),JPKS)
    ENDDO
  ENDDO
ENDDO
CALL GSTATS(1801,1)
CALL GSTATS(140,1)

CALL GSTATS(190,0)
CALL SUTRLE(ZPNMG)
CALL GSTATS(190,1)

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
    F%REPSNM(ICOUNT) = REAL(SQRT(REAL(JN*JN-IM*IM,JPKD)/&
     &REAL(4*JN*JN-1,JPKD)),JPKS)
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
F%RLAPIN(:) = 0.0_JPRB
F%RLAPIN(0) = 0._JPRB
F%RLAPIN(-1) = 0.0_JPRB
DO JN=1,R%NSMAX+2
  F%RLAPIN(JN)=REAL(-(REAL(RA,JPKD)*REAL(RA,JPKD))/REAL(JN*(JN+1),JPKD),JPKS)
ENDDO

DEALLOCATE(ZPNMG)

!     ------------------------------------------------------------------
9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

END SUBROUTINE SULEG
END MODULE SULEG_MOD
