! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FSC_MOD
CONTAINS
SUBROUTINE FSC(KF_UV,KF_SCALARS,KF_SCDERS,&
 &             KST_UV,KST_SC,KST_NSDERS,KST_EWDERS,KST_UVDERS)

!**** *FSC - Division by a*cos(theta), east-west derivatives

!     Purpose.
!     --------
!        In Fourier space divide u and v and all north-south
!        derivatives by a*cos(theta). Also compute east-west derivatives
!        of u,v,thermodynamic, passiv scalar variables and surface
!        pressure.

!**   Interface.
!     ----------
!        CALL FSC(..)
!        Explicit arguments :  PUV     - u and v
!        --------------------  PSCALAR - scalar valued varaibles
!                              PNSDERS - N-S derivative of S.V.V.
!                              PEWDERS - E-W derivative of S.V.V.
!                              PUVDERS - E-W derivative of u and v
!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03 (From SC2FSC)

!     ------------------------------------------------------------------

USE PARKIND_ECTRANS ,ONLY : JPIM     ,JPRBT

USE TPM_TRANS       ,ONLY : LUVDER, LATLON, ZGTF
USE TPM_DISTR       ,ONLY : D, MYSETW,  MYPROC, NPROC, D_NUMP,  D_NPTRLS, D_NSTAGTF
USE TPM_FIELDS      ,ONLY : F
USE TPM_GEOMETRY    ,ONLY : G, G_NMEN
USE TPM_FLT         ,ONLY : S
USE TPM_GEN         ,ONLY : NOUT
USE TPM_DIM         ,ONLY : R
!

IMPLICIT NONE
INTEGER(KIND=JPIM) :: KGL
INTEGER(KIND=JPIM) , INTENT(IN) :: KF_UV,KF_SCALARS,KF_SCDERS
INTEGER(KIND=JPIM) , INTENT(IN) :: KST_UV, KST_SC, KST_NSDERS, KST_EWDERS, KST_UVDERS

REAL(KIND=JPRBT) , POINTER :: PUV(:,:)
REAL(KIND=JPRBT) , POINTER :: PSCALAR(:,:)
REAL(KIND=JPRBT) , POINTER :: PNSDERS(:,:)
REAL(KIND=JPRBT) , POINTER :: PEWDERS(:,:)
REAL(KIND=JPRBT) , POINTER :: PUVDERS(:,:)

REAL(KIND=JPRBT) :: ZACHTE(R%NDGL),ZMUL, ZACHTE2, ZSHIFT, ZPI
REAL(KIND=JPRBT) :: ZAMP, ZPHASE
INTEGER(KIND=JPIM) :: IMEN,ISTAGTF


INTEGER(KIND=JPIM) :: JLON,JF,IGLG,II,IR,JM

INTEGER(KIND=JPIM) :: IBEG,IEND,IINC
!DEBUGGING:
integer :: i,J,maxi,maxj
real :: maxv

!     ------------------------------------------------------------------

IF(MYPROC > NPROC/2)THEN
  IBEG=1
  IEND=D%NDGL_FS
  IINC=1
ELSE
  IBEG=D%NDGL_FS
  IEND=1
  IINC=-1
ENDIF

IF( KF_UV > 0 ) THEN
  PUV     => ZGTF(2*KST_UV-1:2*(KST_UV+2*KF_UV-1),:)
ENDIF
PSCALAR => ZGTF(2*KST_SC-1:2*(KST_SC+KF_SCALARS-1),:)
IF( KF_SCDERS > 0 ) THEN
  PNSDERS => ZGTF(2*KST_nsders-1:2*(KST_nsders+KF_SCDERS-1),:)
  PEWDERS => ZGTF(2*KST_ewders-1:2*(KST_ewders+KF_SCDERS-1),:)
ENDIF
IF (LUVDER) THEN
  PUVDERS => ZGTF(2*KST_uvders-1:2*(KST_uvders+2*KF_UV-1),:)
ENDIF

ZACHTE(:)  = F%RACTHE(:)

#ifdef ACCGPU
!$ACC DATA &
!$ACC& COPYIN(ZACHTE,IBEG,IEND,IINC,KF_SCALARS,KF_UV,KF_SCDERS,MYSETW) &
!$ACC& PRESENT(ZGTF,D_NPTRLS,G_NMEN,D_NSTAGTF) &
!$ACC& PRESENT(PUV,PSCALAR,PNSDERS,PEWDERS,PUVDERS)
#endif
#ifdef OMPGPU
!$OMP TARGET DATA MAP(PRESENT,ALLOC:ZGTF) &
!$OMP& MAP(ALLOC:PUV,PSCALAR,PNSDERS,PEWDERS,PUVDERS)
#endif

DO KGL=IBEG,IEND,IINC

IGLG    = D_NPTRLS(MYSETW)+KGL-1
IMEN    = G_NMEN(IGLG)
ISTAGTF = D_NSTAGTF(KGL)
ZACHTE2  = ZACHTE(IGLG)

!     ------------------------------------------------------------------

!*       1.    DIVIDE U V AND N-S DERIVATIVES BY A*COS(THETA)
!              ----------------------------------------------


!*       1.1      U AND V.

IF(KF_UV > 0) THEN
#ifdef OMPGPU
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) DEFAULT(NONE) SHARED(ISTAGTF,IMEN,KF_UV,PUV,ZACHTE2)
#endif
#ifdef ACCGPU
  !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE) &
  !$ACC& COPYIN(KF_UV,IMEN,ISTAGTF,ZACHTE2) &
  !$ACC& PRESENT(PUV)
#endif
  DO JLON=ISTAGTF+1,ISTAGTF+2*IMEN+1
    DO JF=1,2*KF_UV
      PUV(2*JF-1,JLON) = PUV(2*JF-1,JLON)*ZACHTE2
      PUV(2*JF,  JLON) = PUV(2*JF  ,JLON)*ZACHTE2
    ENDDO
  ENDDO
ENDIF

!*      1.2      N-S DERIVATIVES

IF(KF_SCDERS > 0)THEN
#ifdef OMPGPU
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) DEFAULT(NONE) SHARED(ISTAGTF,IMEN,KF_SCALARS,PNSDERS,ZACHTE2)
#endif
#ifdef ACCGPU
  !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE) &
  !$ACC& COPYIN(KF_SCALARS,IMEN,ISTAGTF,ZACHTE2) &
  !$ACC& PRESENT(PNSDERS)
#endif
  DO JLON=ISTAGTF+1,ISTAGTF+2*IMEN+1
    DO JF=1,KF_SCALARS
      PNSDERS(2*JF-1,JLON) = PNSDERS(2*JF-1,JLON)*ZACHTE2
      PNSDERS(2*JF,  JLON) = PNSDERS(2*JF,  JLON)*ZACHTE2
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

!*       2.    EAST-WEST DERIVATIVES
!              ---------------------

!*       2.1      U AND V.

IF(LUVDER)THEN
#ifdef OMPGPU
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) DEFAULT(NONE) PRIVATE(IR) SHARED(ISTAGTF,IMEN,KF_UV,PUVDERS,ZACHTE2,PUV)
#endif
#ifdef ACCGPU
  !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(IR) DEFAULT(NONE) &
  !$ACC& COPYIN(KF_UV,IMEN,ISTAGTF,ZACHTE2) &
  !$ACC& PRESENT(PUV,PUVDERS)
#endif
  DO JM=0,IMEN
    DO JF=1,2*KF_UV
      IR = ISTAGTF+2*JM+1
      PUVDERS(2*JF-1,IR) = -PUV(2*JF,IR)*ZACHTE2*REAL(JM,JPRBT)
      PUVDERS(2*JF,  IR) =  PUV(2*JF-1,IR)*ZACHTE2*REAL(JM,JPRBT)
    ENDDO
  ENDDO
ENDIF

!*       2.2     SCALAR VARIABLES

IF(KF_SCDERS > 0)THEN
#ifdef OMPGPU
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) DEFAULT(NONE) PRIVATE(IR) SHARED(ISTAGTF,IMEN,KF_SCALARS,PEWDERS,ZACHTE2,PSCALAR)
#endif
#ifdef ACCGPU
  !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(IR) DEFAULT(NONE) &
  !$ACC& COPYIN(KF_SCALARS,IMEN,ISTAGTF,ZACHTE2) &
  !$ACC& PRESENT(PEWDERS,PSCALAR)
#endif
  DO JM=0,IMEN
    DO JF=1,KF_SCALARS
      IR = ISTAGTF+2*JM+1
      PEWDERS(2*JF-1,IR) = -PSCALAR(2*JF,IR)*ZACHTE2*REAL(JM,JPRBT)
      PEWDERS(2*JF,  IR) =  PSCALAR(2*JF-1,IR)*ZACHTE2*REAL(JM,JPRBT)
    ENDDO
  ENDDO
ENDIF

enddo
#ifdef ACCGPU
!$ACC END DATA
#endif
#ifdef OMPGPU
!$OMP END TARGET DATA
#endif
!     ------------------------------------------------------------------

END SUBROUTINE FSC
END MODULE FSC_MOD
