! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE ELTDIR_CTL_MOD
CONTAINS
SUBROUTINE ELTDIR_CTL(KF_FS,KF_UV,KF_SCALARS, &
 & PSPVOR,PSPDIV,PSPSCALAR, &
 & PSPSC3A,PSPSC3B,PSPSC2, &
 & KFLDPTRUV,KFLDPTRSC,PSPMEANU,PSPMEANV,AUX_PROC)

!**** *ELTDIR_CTL* - Control routine for direct Legendre transform

!     Purpose.
!     --------
!        Direct Legendre transform

!**   Interface.
!     ----------
!     CALL ELTDIR_CTL(...)

!     Explicit arguments :
!     --------------------
!     KF_FS      - number of fields in Fourier space
!     KF_UV      - local number of spectral u-v fields
!     KF_SCALARS - local number of scalar spectral fields
!     PSPVOR(:,:) - spectral vorticity (output)
!     PSPDIV(:,:) - spectral divergence (output)
!     PSPSCALAR(:,:) - spectral scalarvalued fields (output)
!     KFLDPTRUV(:) - field pointer for vorticity and divergence (input)
!     KFLDPTRSC(:) - field pointer for scalarvalued fields (input)
!     PSPMEANU(:),PSPMEANV(:) - mean winds

!   R. El Khatib 02-Jun-2022 Optimization/Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_GEN         ,ONLY : LALLOPERM
USE TPM_TRANS       ,ONLY : FOUBUF, FOUBUF_IN
USE TPM_DISTR       ,ONLY : D

USE ELTDIR_MOD      ,ONLY : ELTDIR
USE EUVTVD_COMM_MOD , ONLY : EUVTVD_COMM
USE TRLTOM_MOD      ,ONLY : TRLTOM
USE MPL_MODULE

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KF_FS,KF_UV,KF_SCALARS
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPSC3A(:,:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPSC3B(:,:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPSC2(:,:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KFLDPTRUV(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KFLDPTRSC(:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPMEANU(:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(OUT) :: PSPMEANV(:)
EXTERNAL AUX_PROC
OPTIONAL AUX_PROC

INTEGER(KIND=JPIM) :: JM,IM,IBLEN,ILED2,INUL
REAL(KIND=JPRB) :: ZDUM
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------

! Transposition from Fourier space distribution to spectral space distribution

IF (LHOOK) CALL DR_HOOK('ELTDIR_CTL_MOD:ELTDIR_CTL',0,ZHOOK_HANDLE)
IBLEN = D%NLENGT0B*2*KF_FS
IF (ALLOCATED(FOUBUF)) THEN
  IF (MAX(1,IBLEN) > SIZE(FOUBUF)) THEN
    DEALLOCATE(FOUBUF)
    ALLOCATE(FOUBUF(MAX(1,IBLEN)))
  ENDIF
ELSE
  ALLOCATE(FOUBUF(MAX(1,IBLEN)))
  FOUBUF(1)=0._JPRB ! enforce allocation here
ENDIF
CALL GSTATS(153,0)
CALL TRLTOM(FOUBUF_IN,FOUBUF,2*KF_FS)
CALL GSTATS(153,1)
IF (.NOT.LALLOPERM) DEALLOCATE(FOUBUF_IN)

! Periodization of auxiliary fields in y direction

IF (PRESENT(AUX_PROC)) THEN
  CALL AUX_PROC(ZDUM,FOUBUF,2*KF_FS,1,IBLEN,0,D%NUMP,.FALSE.,&
   & INUL,D%NPROCL,D%NSTAGT0B,D%NPNTGTB1)
ENDIF

! Direct Legendre transform

ILED2 = 2*KF_FS
CALL GSTATS(1645,0)
IF (KF_FS>0) THEN
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JM,IM)
  DO JM=1,D%NUMP
    IM = D%MYMS(JM)
    CALL ELTDIR(IM,JM,KF_FS,KF_UV,KF_SCALARS,ILED2, &
     & PSPVOR,PSPDIV,PSPSCALAR,&
     & PSPSC3A,PSPSC3B,PSPSC2 , &
     & KFLDPTRUV,KFLDPTRSC,PSPMEANU,PSPMEANV)
  ENDDO
!$OMP END PARALLEL DO
  IF (KF_UV > 0) THEN
    CALL EUVTVD_COMM(KF_UV,PSPMEANU,PSPMEANV,KFLDPTRUV)
  ENDIF
ENDIF
CALL GSTATS(1645,1)
  
IF (.NOT.LALLOPERM) DEALLOCATE(FOUBUF)
IF (LHOOK) CALL DR_HOOK('ELTDIR_CTL_MOD:ELTDIR_CTL',1,ZHOOK_HANDLE)

!     -----------------------------------------------------------------

END SUBROUTINE ELTDIR_CTL
END MODULE ELTDIR_CTL_MOD
