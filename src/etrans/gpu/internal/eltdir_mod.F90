#define ALIGN(I, A) (((I)+(A)-1)/(A)*(A))
MODULE ELTDIR_MOD
  USE BUFFERED_ALLOCATOR_MOD

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: ELTDIR, ELTDIR_HANDLE, PREPARE_ELTDIR

  TYPE ELTDIR_HANDLE
    TYPE(ALLOCATION_RESERVATION_HANDLE) :: HFFT_AND_VODI
  END TYPE

CONTAINS
  FUNCTION PREPARE_ELTDIR(ALLOCATOR,KF_FS,KF_UV) RESULT(HELTDIR)
    USE PARKIND_ECTRANS, ONLY: JPIM, JPRBT, JPRD
    USE TPM_DISTR, ONLY: D
    USE TPM_DIM, ONLY: R
    USE TPMALD_DIM      ,ONLY : RALD
    USE ISO_C_BINDING
    USE LEINV_MOD

    IMPLICIT NONE

    TYPE(BUFFERED_ALLOCATOR), INTENT(INOUT) :: ALLOCATOR
    INTEGER(KIND=JPIM), INTENT(IN) ::KF_FS, KF_UV

    TYPE(ELTDIR_HANDLE) :: HELTDIR

    INTEGER(KIND=C_SIZE_T) :: IALLOC_SZ
    REAL(KIND=JPRBT) :: ZPRBT_DUMMY
    
    ! ZFFT
    IALLOC_SZ = ALIGN((RALD%NDGLSUR+R%NNOEXTZG)*D%NUMP*2*KF_FS*SIZEOF(ZPRBT_DUMMY), 128)
    ! ZVODI
    IALLOC_SZ = IALLOC_SZ+ALIGN((RALD%NDGLSUR+R%NNOEXTZG)*D%NUMP*MAX(4*KF_UV,1)*SIZEOF(ZPRBT_DUMMY), 128)
    HELTDIR%HFFT_AND_VODI = RESERVE(ALLOCATOR, IALLOC_SZ)
    
  END FUNCTION PREPARE_ELTDIR

SUBROUTINE ELTDIR(ALLOCATOR,HELTDIR,KF_FS,KF_UV,KF_SCALARS,FOUBUF,&
 & PSPVOR,PSPDIV,PSPSCALAR,PSPSC3A,PSPSC3B,PSPSC2,&
 & KFLDPTRUV,KFLDPTRSC,PSPMEANU,PSPMEANV)


USE ISO_C_BINDING

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DIM         ,ONLY : R
USE TPM_DISTR       ,ONLY : D, d_myms, d_nump
USE TPMALD_DIM      ,ONLY : RALD

USE EPRFI2B_MOD      ,ONLY : EPRFI2B
USE ELEDIR_MOD      ,ONLY : ELEDIR
USE EUVTVD_MOD
USE EUVTVD_COMM_MOD
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

TYPE(BUFFERED_ALLOCATOR), INTENT(IN) :: ALLOCATOR
TYPE(ELTDIR_HANDLE), INTENT(IN) :: HELTDIR

INTEGER(KIND=JPIM),INTENT(IN)   :: KF_FS,KF_UV,KF_SCALARS
REAL(KIND=JPRB), INTENT(IN) :: FOUBUF(:)
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
INTEGER(KIND=JPIM) :: IFLD, IR, J, IM, JM
INTEGER(KIND=JPIM) :: IUS,IVS,IVORS,IDIVS, IUE, IVE, IVORE, IDIVE

REAL(KIND=JPRB), POINTER :: ZFFT(:,:,:), ZFFT_L(:)
REAL(KIND=JPRB), POINTER :: ZVODI(:,:,:), ZVODI_L(:)
INTEGER(KIND=C_SIZE_T) :: IALLOC_SZ, IALLOC_POS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ELTDIR_MOD:ELTDIR',0,ZHOOK_HANDLE)

! ZFFT(RALD%NDGLSUR+R%NNOEXTZG,2*KF_FS,D%NUMP)
IALLOC_POS = 1
IALLOC_SZ = ALIGN((RALD%NDGLSUR+R%NNOEXTZG)*D%NUMP*2*KF_FS*SIZEOF(ZFFT_L(1)), 128)
CALL ASSIGN_PTR(ZFFT_L, GET_ALLOCATION(ALLOCATOR, HELTDIR%HFFT_AND_VODI),&
    & IALLOC_POS, IALLOC_SZ)
CALL C_F_POINTER(C_LOC(ZFFT_L), ZFFT, (/ RALD%NDGLSUR+R%NNOEXTZG,D%NUMP,2*KF_FS /))
IALLOC_POS = IALLOC_POS + IALLOC_SZ
! ZVODI(RALD%NDGLSUR+R%NNOEXTZG,MAX(4*KF_UV,1),D%NUMP)
IALLOC_SZ = ALIGN((RALD%NDGLSUR+R%NNOEXTZG)*D%NUMP*MAX(4*KF_UV,1)*SIZEOF(ZVODI_L(1)), 128)
CALL ASSIGN_PTR(ZVODI_L, GET_ALLOCATION(ALLOCATOR, HELTDIR%HFFT_AND_VODI),&
    & IALLOC_POS, IALLOC_SZ)
CALL C_F_POINTER(C_LOC(ZVODI_L), ZVODI, (/ RALD%NDGLSUR+R%NNOEXTZG,D%NUMP,MAX(4*KF_UV,1) /))
IALLOC_POS = IALLOC_POS + IALLOC_SZ

!*     1.    PREPARE WORK ARRAYS.
!            --------------------

CALL EPRFI2B(KF_FS,ZFFT,FOUBUF)

!*     2.    PERIODICIZATION IN Y DIRECTION
!            ------------------------------

IF(R%NNOEXTZG>0) THEN
  CALL ABORT('NNOEXTZG>0 not supported on GPU')
ENDIF

!*     3.    DIRECT LEGENDRE TRANSFORM.
!            --------------------------

CALL ELEDIR(ALLOCATOR,ZFFT)

!*     4.    COMPUTE VORTICITY AND DIVERGENCE AND STORE MEAN WIND ON TASK OWNING WAVE 0
!            --------------------------------------------------------------------------

#ifdef ACCGPU
    !$ACC DATA COPYOUT(PSPVOR,PSPDIV) IF(KF_UV > 0)
    !$ACC DATA COPYOUT(PSPSCALAR) IF(PRESENT(PSPSCALAR) .AND. KF_SCALARS > 0)
    !$ACC DATA COPYOUT(PSPSC2) IF(PRESENT(PSPSC2))
    !$ACC DATA COPYOUT(PSPSC3A) IF((PRESENT(PSPSC3A)))
    !$ACC DATA COPYOUT(PSPSC3B) IF((PRESENT(PSPSC3B)))
#endif

IF( KF_UV > 0 ) THEN
  IUS = 1
  IUE = 2*KF_UV
  IVS = 2*KF_UV+1
  IVE = 4*KF_UV
  IVORS = 1
  IVORE = 2*KF_UV
  IDIVS = 2*KF_UV+1
  IDIVE = 4*KF_UV
  CALL EUVTVD(KF_UV,ZFFT(:,:,IUS:IUE),ZFFT(:,:,IVS:IVE),&
    & ZVODI(:,:,IVORS:IVORE),ZVODI(:,:,IDIVS:IDIVE))
    
  DO JM=1,D%NUMP
    IM = D%MYMS(JM)

    CALL EUVTVD_COMM(IM,JM,KF_UV,KFLDPTRUV,ZFFT(:,:,IUS:IUE), &
     & ZFFT(:,:,IVS:IVE), &
     & PSPMEANU,PSPMEANV)

  ENDDO
  
ENDIF

!*       5.    UPDATE SPECTRAL ARRAYS.
!              -----------------------

CALL EUPDSP(KF_UV,KF_SCALARS,ZFFT,ZVODI, &
 & PSPVOR,PSPDIV,PSPSCALAR,PSPSC3A,PSPSC3B,PSPSC2,KFLDPTRUV,KFLDPTRSC)

#ifdef ACCGPU
    !$ACC END DATA
    !$ACC END DATA
    !$ACC END DATA
    !$ACC END DATA
    !$ACC END DATA
#endif


IF (LHOOK) CALL DR_HOOK('ELTDIR_MOD:ELTDIR',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE ELTDIR
END MODULE ELTDIR_MOD
