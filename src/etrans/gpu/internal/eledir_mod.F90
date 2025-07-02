MODULE ELEDIR_MOD
CONTAINS
SUBROUTINE ELEDIR(ALLOCATOR,PFFT,PFFT_OUT)

!**** *LEINV* - Inverse Legendre transform.

!     Purpose.
!     --------
!        Inverse Legendre tranform of all variables(kernel).

!**   Interface.
!     ----------
!        CALL LEINV(...)

!        Explicit arguments :  KM - zonal wavenumber (input-c)
!        --------------------  KFC - number of fields to tranform (input-c)
!                              PIA - spectral fields
!                              for zonal wavenumber KM (input)
!                              PAOA1 - antisymmetric part of Fourier
!                              fields for zonal wavenumber KM (output)
!                              PSOA1 - symmetric part of Fourier
!                              fields for zonal wavenumber KM (output)
!                              PLEPO - Legendre polonomials for zonal
!                              wavenumber KM (input-c)

!        Implicit arguments :  None.
!        --------------------

!     Method.
!     -------

!     Externals.   MXMAOP - calls SGEMVX (matrix multiply)
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From LEINV in IFS CY22R1
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        R. El Khatib 01-Sep-2015 support for FFTW transforms

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM, JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_GEN,                ONLY: NCUR_RESOL
USE TPM_DISTR       ,ONLY : D
USE TPMALD_DIM      ,ONLY : RALD

USE TPM_HICFFT      ,ONLY : EXECUTE_DIR_FFT

USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS

USE ISO_C_BINDING
USE BUFFERED_ALLOCATOR_MOD

!

IMPLICIT NONE

TYPE(BUFFERED_ALLOCATOR), INTENT(IN) :: ALLOCATOR
REAL(KIND=JPRB),  POINTER,  INTENT(INOUT)  :: PFFT(:,:,:)
REAL(KIND=JPRB),  POINTER,  INTENT(INOUT)  :: PFFT_OUT(:,:,:)

INTEGER(KIND=JPIM) :: JLOT
TYPE(C_PTR) :: IPLAN_C2R
REAL (KIND=JPRB), POINTER :: ZFFT_L(:), ZFFT_L_OUT(:)  ! 1D copy
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------

IF (LHOOK) CALL DR_HOOK('ELEDIR_MOD:ELEDIR',0,ZHOOK_HANDLE)

JLOT=UBOUND(PFFT,2)*UBOUND (PFFT,3)

CALL C_F_POINTER(C_LOC(PFFT), ZFFT_L, (/UBOUND(PFFT,1)*UBOUND(PFFT,2)*UBOUND(PFFT,3)/) )
CALL C_F_POINTER(C_LOC(PFFT_OUT), ZFFT_L_OUT, (/UBOUND(PFFT_OUT,1)*UBOUND(PFFT_OUT,2)*UBOUND(PFFT_OUT,3)/) )

IF (JLOT==0) THEN
  IF (LHOOK) CALL DR_HOOK('ELEDIR_MOD:ELEDIR',1,ZHOOK_HANDLE)
  RETURN
ENDIF

!$ACC DATA PRESENT(ZFFT_L,ZFFT_L_OUT,RALD%NLOENS_LAT,RALD%NOFFSETS_LAT)
CALL EXECUTE_DIR_FFT(ZFFT_L(:),ZFFT_L_OUT(:),NCUR_RESOL,-JLOT, &    ! -JLOT to have hicfft make distinction between zonal and meridional direction. Don't worry, abs(JLOT) is used internally ...
    & LOENS=RALD%NLOENS_LAT, &
    & OFFSETS=RALD%NOFFSETS_LAT,ALLOC=ALLOCATOR%PTR)
!$ACC END DATA

IF (LHOOK) CALL DR_HOOK('ELEDIR_MOD:ELEDIR',1,ZHOOK_HANDLE)

END SUBROUTINE ELEDIR
END MODULE ELEDIR_MOD