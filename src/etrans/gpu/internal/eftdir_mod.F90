MODULE EFTDIR_MOD
CONTAINS
SUBROUTINE EFTDIR(ALLOCATOR,PREEL,KF_FS,AUX_PROC)

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
INTEGER(KIND=JPIM), INTENT(IN)  :: KF_FS
REAL(KIND=JPRB),    INTENT(INOUT)  :: PREEL(:)   ! (IRLEN+2)*NDGLG*KF_FS
EXTERNAL AUX_PROC
OPTIONAL AUX_PROC

INTEGER(KIND=JPIM) :: JLOT, IRLEN, JJ
! INTEGER(KIND=JPIB), SAVE :: OFFSETS(2)
! INTEGER(KIND=JPIM), SAVE :: LOENS(1)
integer :: istat
character(len=32) :: cfrmt
REAL(KIND=JPRB) :: ZDUM
INTEGER(KIND=JPIM) :: INUL

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------

IF (LHOOK) CALL DR_HOOK('ELEINV_MOD:ELEINV',0,ZHOOK_HANDLE)

IRLEN=R%NDLON+R%NNOEXTZG


! Periodization of auxiliary fields in x direction
IF(R%NNOEXTZL>0) THEN
  !!! FIXME !!! CALL EXTPER(PREEL,R%NDLON+R%NNOEXTZL,1,R%NDLON,KF_FS,D%NDGL_FS,D%NSTAGTF,0)
  CALL ABORT('FIXME')
ENDIF
IF (PRESENT(AUX_PROC)) THEN
  !!! FIXME !!! CALL AUX_PROC(PREEL,ZDUM,KF_FS,D%NLENGTF,1,D%NDGL_FS,0,.TRUE.,&
  !!! FIXME !!!  & D%NSTAGTF,INUL,INUL,INUL)
  CALL ABORT('FIXME')
ENDIF


! LOENS(1)=IRLEN
JLOT=SIZE(PREEL)/(IRLEN+2)
! compute offsets; TODO: avoid recomputing/putting on device every time.
! DO JJ=1,SIZE(OFFSETS)
  ! OFFSETS(JJ)=(JJ-1)*(IRLEN+2)
! ENDDO

IF (JLOT==0) THEN
  IF (LHOOK) CALL DR_HOOK('ELEINV_MOD:ELEINV',1,ZHOOK_HANDLE)
  RETURN
ENDIF


#ifdef gnarls
! ! debugging
! !$acc data present(preel)
! !$acc update host(preel)
! write (6,*) 'before direct zonal transform : '
! write (6,*) 'shape(preel) = ',shape(preel)
! write (6,*) 'preel = ',preel
! call flush(6)
! !$acc end data
#endif


#ifdef gnarls

! fake fft: only take mean value, on cpu
!$acc data present(preel)
!$acc update host(preel)
DO JJ=1,JLOT
  preel((JJ-1)*(irlen+2)+1)=sum(preel((JJ-1)*(irlen+2)+1:jj*(irlen+2)-2))   ! mean value
  preel((JJ-1)*(irlen+2)+2:jj*(irlen+2))=0.
ENDDO
!$acc update device(preel)
!$acc end data

#else

!write (6,*) __FILE__, __LINE__; call flush(6)
!$ACC DATA PRESENT(PREEL,RALD%NLOENS_LON,RALD%NOFFSETS_LON)
!write (6,*) __FILE__, __LINE__; call flush(6)
CALL EXECUTE_DIR_FFT(PREEL(:),PREEL(:),NCUR_RESOL,JLOT, &
    & LOENS=RALD%NLOENS_LON, &
    & OFFSETS=RALD%NOFFSETS_LON,ALLOC=ALLOCATOR%PTR)
!write (6,*) __FILE__, __LINE__; call flush(6)
!$ACC END DATA
!write (6,*) __FILE__, __LINE__; call flush(6)


#endif


#ifdef gnarls
! ! debugging
! !$acc data present(preel)
! !$acc update host(preel)
! write (6,*) 'before direct zonal transform : '
! write (6,*) 'shape(preel) = ',shape(preel)
! write (6,*) 'preel = ',preel
! call flush(6)
! !$acc end data
#endif

IF (LHOOK) CALL DR_HOOK('EFTDIR_MOD:EFTDIR',1,ZHOOK_HANDLE)

END SUBROUTINE EFTDIR
END MODULE EFTDIR_MOD