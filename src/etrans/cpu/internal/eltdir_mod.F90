MODULE ELTDIR_MOD
CONTAINS
SUBROUTINE ELTDIR(KM,KMLOC,KF_FS,KF_UV,KF_SCALARS,KLED2,&
 & PSPVOR,PSPDIV,PSPSCALAR,PSPSC3A,PSPSC3B,PSPSC2,&
 & KFLDPTRUV,KFLDPTRSC,PSPMEANU,PSPMEANV)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DIM         ,ONLY : R
USE TPM_DISTR       ,ONLY : D
USE TPMALD_DIM      ,ONLY : RALD

USE EPRFI2_MOD      ,ONLY : EPRFI2
USE ELEDIR_MOD      ,ONLY : ELEDIR
USE EUVTVD_MOD
USE EUPDSP_MOD      ,ONLY : EUPDSP
USE EXTPER_MOD      ,ONLY : EXTPER

!
!**** *ELTDIR* - Control of Direct Legendre transform step

!     Purpose.
!     --------
!        Tranform from Fourier space to spectral space, compute
!        vorticity and divergence.

!**   Interface.
!     ----------
!        *CALL* *ELTDIR(...)*

!        Explicit arguments :
!        --------------------  KM     - zonal wavenumber
!                              KMLOC  - local zonal wavenumber

!        Implicit arguments :  None
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------
!         EPRFI2      - prepares the Fourier work arrays for model variables
!         ELEDIR      - direct Legendre transform
!         EUVTVD      -
!         EUPDSP      - updating of spectral arrays (fields)
!         EUVTVD_COMM -
!         EXTPER      -


!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 87-11-24
!        Modified : 91-07-01 Philippe Courtier/Mats Hamrud - Rewrite
!                            for uv formulation
!        Modified 93-03-19 D. Giard - CDCONF='T' for tendencies
!        Modified 93-11-18 M. Hamrud - use only one Fourier buffer
!        Modified 94-04-06 R. El khatib Full-POS implementation
!        M.Hamrud  : 94-11-01 New conf 'G' - vor,div->vor,div
!                             instead of u,v->vor,div
!        MPP Group : 95-10-01 Support for Distributed Memory version
!        K. YESSAD (AUGUST 1996):
!               - Legendre transforms for transmission coefficients.
!        Modified : 04/06/99 D.Salmond : change order of AIA and SIA
!            01-03-14 G. Radnoti aladin version
!     01-08-27 : R. El Khatib Fix for NPROMATR /= 0
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        F. Vana + NEC 28-Apr-2009 MPI-OpenMP fix
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        N. Lopes & R. El Khatib 15-Jun-2012 Scalability enhancement
!        R. El Khatib 02-Jun-2022 Optimization/Cleaning
!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)   :: KM
INTEGER(KIND=JPIM),INTENT(IN)   :: KMLOC
INTEGER(KIND=JPIM),INTENT(IN)   :: KF_FS,KF_UV,KF_SCALARS,KLED2

REAL(KIND=JPRB)  ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
REAL(KIND=JPRB)  ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
REAL(KIND=JPRB)  ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPSC2(:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPSC3A(:,:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPSC3B(:,:,:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KFLDPTRUV(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KFLDPTRSC(:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPMEANU(:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPMEANV(:)

INTEGER(KIND=JPIM) :: IFC, IINDEX(2*KF_FS), JF, JDIM
INTEGER(KIND=JPIM) :: IFLD, IR, J
INTEGER(KIND=JPIM) :: IUS,IVS,IVORS,IDIVS

REAL(KIND=JPRB) :: ZFFT(RALD%NDGLSUR+R%NNOEXTZG,KLED2,D%NUMP)
REAL(KIND=JPRB) :: ZVODI(RALD%NDGLSUR+R%NNOEXTZG,MAX(4*KF_UV,1),D%NUMP)

! Only if R%NNOEXTZG > 0 :
REAL(KIND=JPRB) :: ZFFT2(KLED2,(RALD%NDGLSUR+R%NNOEXTZG)*MIN(1,R%NNOEXTZG))

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ELTDIR_MOD:ELTDIR',0,ZHOOK_HANDLE)

IUS = 1
IVS = 2*KF_UV+1
IVORS = IUS
IDIVS = IVS
IFC = 2*KF_FS

!*     1.    PREPARE WORK ARRAYS.
!            --------------------

CALL EPRFI2(KM,KMLOC,KF_FS,ZFFT(:,:,KMLOC))

!*     2.    PERIODICIZATION IN Y DIRECTION
!            ------------------------------

IF(R%NNOEXTZG>0) THEN
  DO JF = 1,IFC
    DO JDIM = 1,R%NDGL
      ZFFT2(JF,JDIM)=ZFFT(JDIM,JF,KMLOC)
    ENDDO
  ENDDO
  IINDEX(1)=0
  CALL EXTPER(ZFFT2(:,:),R%NDGL+R%NNOEXTZG,1,R%NDGL,IFC,1,IINDEX,0)
  DO JF = 1,IFC
    DO JDIM = 1,R%NDGL+R%NNOEXTZG
      ZFFT(JDIM,JF,KMLOC) = ZFFT2(JF,JDIM)
    ENDDO
  ENDDO
ENDIF

!*     3.    DIRECT LEGENDRE TRANSFORM.
!            --------------------------

CALL ELEDIR(KM,IFC,KLED2,ZFFT(:,:,KMLOC))

!*     4.    COMPUTE VORTICITY AND DIVERGENCE AND STORE MEAN WIND ON TASK OWNING WAVE 0
!            --------------------------------------------------------------------------

IF( KF_UV > 0 ) THEN
  CALL EUVTVD(KM,KMLOC,KF_UV,ZFFT(:,IUS:,KMLOC),ZFFT(:,IVS:,KMLOC),&
   & ZVODI(:,IVORS:,KMLOC),ZVODI(:,IDIVS:,KMLOC))
  IF (KM == 0) THEN
    IF (PRESENT(KFLDPTRUV)) THEN
      DO J = 1, KF_UV
        IR = 2*J-1
        IFLD=KFLDPTRUV(J)
        PSPMEANU(IFLD)=ZFFT(1,IUS-1+IR,KMLOC)
        PSPMEANV(IFLD)=ZFFT(1,IVS-1+IR,KMLOC)
      ENDDO
    ELSE
      DO J = 1, KF_UV
        IR = 2*J-1
        PSPMEANU(J)=ZFFT(1,IUS-1+IR,KMLOC)
        PSPMEANV(J)=ZFFT(1,IVS-1+IR,KMLOC)
      ENDDO
    ENDIF
  ENDIF
ENDIF

!*       5.    UPDATE SPECTRAL ARRAYS.
!              -----------------------

CALL EUPDSP(KM,KF_UV,KF_SCALARS,ZFFT(:,:,KMLOC),ZVODI(:,:,KMLOC), &
 & PSPVOR,PSPDIV,PSPSCALAR,PSPSC3A,PSPSC3B,PSPSC2,KFLDPTRUV,KFLDPTRSC)

IF (LHOOK) CALL DR_HOOK('ELTDIR_MOD:ELTDIR',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE ELTDIR
END MODULE ELTDIR_MOD
