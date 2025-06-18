MODULE ELEDIR_MOD
CONTAINS
SUBROUTINE ELEDIR(ALLOCATOR,PFFT)

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

USE PARKIND1  ,ONLY : JPIM, JPRB, JPIB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_GEN,                ONLY: NCUR_RESOL
USE TPM_DISTR       ,ONLY : D
USE TPM_DIM         ,ONLY : R
USE TPMALD_DIM      ,ONLY : RALD
USE TPMALD_FFT      ,ONLY : TALD

USE TPM_HICFFT      ,ONLY : EXECUTE_DIR_FFT

USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS

USE ISO_C_BINDING
USE BUFFERED_ALLOCATOR_MOD

!

IMPLICIT NONE

TYPE(BUFFERED_ALLOCATOR), INTENT(IN) :: ALLOCATOR
REAL(KIND=JPRB),    INTENT(INOUT)  :: PFFT(:,:,:)

INTEGER(KIND=JPIM) :: IRLEN, ICLEN, JLOT, JJ
!INTEGER(KIND=JPIM) :: IPLAN_C2R
TYPE(C_PTR) :: IPLAN_C2R
REAL (KIND=JPRB)   :: ZSCAL
REAL (KIND=JPRB), POINTER :: ZFFT_L(:)  ! 1D copy
INTEGER(KIND=JPIB) :: OFFSETS(2)
INTEGER(KIND=JPIM) :: LOENS(1)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------

IF (LHOOK) CALL DR_HOOK('ELEDIR_MOD:ELEDIR',0,ZHOOK_HANDLE)

IRLEN=R%NDGL+R%NNOEXTZG
LOENS(1)=IRLEN
JLOT=UBOUND(PFFT,2)*UBOUND (PFFT,3)

! compute offsets; TODO: avoid recomputing/putting on device every time.
DO JJ=1,SIZE(OFFSETS)
  OFFSETS(JJ)=(JJ-1)*(IRLEN+2)
ENDDO

write (6,*) __FILE__,__LINE__; call flush(6)
write (6,*)'  JLOT = ',JLOT; call flush(6)
write (6,*)'  shape(PFFT) = ',shape(PFFT); call flush(6)

CALL C_F_POINTER(C_LOC(PFFT), ZFFT_L, (/UBOUND(PFFT,1)*UBOUND(PFFT,2)*UBOUND(PFFT,3)/) )

#ifndef gnarls
! debugging
!$acc data present(zfft_l)
!$acc update host(zfft_l)
write (6,*) 'before meridional transform:'
write (6,*) 'zfft_l = '
write (6,'(8F10.3)') zfft_l
call flush(6)
!$acc end data
#endif

IF (JLOT==0) THEN
  IF (LHOOK) CALL DR_HOOK('ELEDIR_MOD:ELEDIR',1,ZHOOK_HANDLE)
  RETURN
ENDIF

#ifdef gnarls

! fake fft: only take mean value, on cpu
!$acc data present(zfft_l)
!$acc update host(zfft_l)
DO JJ=1,JLOT
  zfft_l((JJ-1)*(irlen+2)+1)=sum(zfft_l((JJ-1)*(irlen+2)+1:jj*(irlen+2)-2))   ! mean value
  zfft_l((JJ-1)*(irlen+2)+2:jj*(irlen+2))=0.
ENDDO
!$acc update device(zfft_l)
!$acc end data

#else

!$ACC DATA PRESENT(PFFT) COPYIN(LOENS,OFFSETS)
CALL EXECUTE_DIR_FFT(ZFFT_L(:),ZFFT_L(:),NCUR_RESOL,-JLOT, &    ! -JLOT to have hicfft make distinction between zonal and meridional direction. Don't worry, abs(JLOT) is used internally ...
    & LOENS=LOENS, &
    & OFFSETS=OFFSETS,ALLOC=ALLOCATOR%PTR)
!$ACC END DATA

#endif


#ifndef gnarls
! debugging
!$acc data present(zfft_l)
!$acc update host(zfft_l)
write (6,*) 'after meridional transform:'
write (6,*) 'zfft_l = '
write (6,'(8F10.3)') zfft_l
call flush(6)
!$acc end data
#endif

IF (LHOOK) CALL DR_HOOK('ELEDIR_MOD:ELEDIR',1,ZHOOK_HANDLE)

END SUBROUTINE ELEDIR
END MODULE ELEDIR_MOD