MODULE EFSC_MOD
CONTAINS
SUBROUTINE EFSC(PREEL,KF_UV,KF_SCALARS,KF_SCDERS,KF_FS)
!SUBROUTINE EFSC(KF_UV,KF_SCALARS,KF_SCDERS,&
! & PUV,PSCALAR,PNSDERS,PEWDERS,PUVDERS)

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
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_TRANS       ,ONLY : LUVDER, LVORGP, LDIVGP
USE TPM_DISTR       ,ONLY : D, MYSETW, D_NPTRLS, D_NSTAGTF
USE TPM_GEOMETRY    ,ONLY : G, G_NMEN, G_NLOEN_MAX
USE TPMALD_GEO      ,ONLY : GALD
!

IMPLICIT NONE

REAL(KIND=JPRB) , INTENT(INOUT) :: PREEL(:)
INTEGER(KIND=JPIM) , INTENT(IN) :: KF_UV,KF_SCALARS,KF_SCDERS, KF_FS

INTEGER(KIND=JPIM) :: IMEN,ISTAGTF

INTEGER(KIND=JPIM) :: JF,IGLG,II,IR,JM,JGL
REAL(KIND=JPRB) :: ZIM
INTEGER(KIND=JPIM) :: I_UV_OFFSET, I_SC_OFFSET, I_SCDERS_OFFSET, I_UVDERS_OFFSET, IST
INTEGER(KIND=JPIM) :: IOFF_LAT,IOFF_UV,IOFF_UV_EWDER, IOFF_SCALARS, IOFF_SCALARS_EWDER
REAL(KIND=JPRB)    :: RET_REAL, RET_COMPLEX
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EFSC_MOD:EFSC',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------


IF(KF_UV > 0 .OR. KF_SCDERS > 0) THEN
  IST = 0
  IF(LVORGP) THEN
    IST = IST+KF_UV
  ENDIF
  IF(LDIVGP) THEN
    IST = IST+KF_UV
  ENDIF
  I_UV_OFFSET=IST

  IST = IST+2*KF_UV
  I_SC_OFFSET=IST

  IST = IST+KF_SCALARS
  !I_NSDERS_OFFSET=IST
  
  IST = IST+KF_SCDERS
  IF(LUVDER) THEN
    I_UVDERS_OFFSET=IST
    IST = IST+2*KF_UV
  ENDIF

  IF(KF_SCDERS > 0) THEN
    I_SCDERS_OFFSET=IST
  ENDIF
ENDIF

#ifdef ACCGPU
!$ACC DATA &
!$ACC& PRESENT(D_NPTRLS,D_NSTAGTF,PREEL,G_NMEN, G_NLOEN_MAX, D)
#endif


!     ------------------------------------------------------------------

!*       2.    EAST-WEST DERIVATIVES
!              ---------------------

!*       2.1      U AND V.

IF (LUVDER) THEN
#ifdef OMPGPU
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO DEFAULT(NONE) SHARED(KF_UVPREEL)
#endif
#ifdef ACCGPU
  !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(IGLG,IOFF_LAT,IOFF_UV,IOFF_UV_EWDER,ZIM,RET_REAL,RET_COMPLEX,JM,JF,JGL) &
  !$ACC& FIRSTPRIVATE(KF_UV,I_UVDERS_OFFSET,I_UV_OFFSET,KF_FS)
#endif
  DO JGL=1,D%NDGL_FS
    DO JF=1,2*KF_UV
      DO JM=0,G_NLOEN_MAX/2
        IGLG = D_NPTRLS(MYSETW)+JGL-1
        IOFF_LAT = KF_FS*D_NSTAGTF(JGL)
        IOFF_UV = IOFF_LAT+(I_UV_OFFSET+JF-1)*(D_NSTAGTF(JGL+1)-D_NSTAGTF(JGL))
        IOFF_UV_EWDER = IOFF_LAT+(I_UVDERS_OFFSET+JF-1)*(D_NSTAGTF(JGL+1)-D_NSTAGTF(JGL))

        RET_REAL = 0.0_JPRBT
        RET_COMPLEX = 0.0_JPRBT

        IF (JM <= G_NMEN(IGLG)) THEN
          ZIM  = REAL(JM,JPRB)*GALD%EXWN

          RET_REAL = &
              & -PREEL(IOFF_UV+2*JM+2)*ZIM
          RET_COMPLEX =  &
              &  PREEL(IOFF_UV+2*JM+1)*ZIM
        ENDIF
        PREEL(IOFF_UV_EWDER+2*JM+1) = RET_REAL
        PREEL(IOFF_UV_EWDER+2*JM+2) = RET_COMPLEX
      ENDDO
    ENDDO
  ENDDO
ENDIF

!*       2.2     SCALAR VARIABLES

IF (KF_SCDERS > 0) THEN
#ifdef OMPGPU
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO DEFAULT(NONE) SHARED(KF_SCALARS,PEWDERS,PSCALAR)
#endif
#ifdef ACCGPU
  !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(IGLG,IOFF_LAT,IOFF_SCALARS_EWDER,IOFF_SCALARS,ZIM,RET_REAL,RET_COMPLEX) &
  !$ACC& FIRSTPRIVATE(KF_SCALARS,I_SCDERS_OFFSET,I_SC_OFFSET,KF_FS)
#endif
  DO JGL=1,D%NDGL_FS
    DO JF=1,KF_SCALARS
      DO JM=0,G_NLOEN_MAX/2
        IGLG = D_NPTRLS(MYSETW)+JGL-1
        IOFF_LAT = KF_FS*D_NSTAGTF(JGL)
        IOFF_SCALARS_EWDER = IOFF_LAT+(I_SCDERS_OFFSET+JF-1)*(D_NSTAGTF(JGL+1)-D_NSTAGTF(JGL))
        IOFF_SCALARS = IOFF_LAT+(I_SC_OFFSET+JF-1)*(D_NSTAGTF(JGL+1)-D_NSTAGTF(JGL))

        RET_REAL = 0.0_JPRBT
        RET_COMPLEX = 0.0_JPRBT

        IF (JM <= G_NMEN(IGLG)) THEN
          ZIM  = REAL(JM,JPRB)*GALD%EXWN

          RET_REAL = &
              & -PREEL(IOFF_SCALARS+2*JM+2)*ZIM
          RET_COMPLEX = &
              &  PREEL(IOFF_SCALARS+2*JM+1)*ZIM
        ENDIF
        ! The rest from G_NMEN(IGLG+1)...MAX is zero truncated
        PREEL(IOFF_SCALARS_EWDER+2*JM+1) = RET_REAL
        PREEL(IOFF_SCALARS_EWDER+2*JM+2) = RET_COMPLEX
      ENDDO
    ENDDO
  ENDDO
ENDIF

#ifdef ACCGPU
!$ACC END DATA
#endif

IF (LHOOK) CALL DR_HOOK('EFSC_MOD:EFSC',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE EFSC
END MODULE EFSC_MOD