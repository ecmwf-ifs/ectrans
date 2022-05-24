! (C) Copyright 2000- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LEINV_MOD
CONTAINS
SUBROUTINE LEINV(KFC,KF_OUT_LT,PIA,PAOA1,PSOA1)

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

!        Implicit arguments :  None.
!        --------------------

!     Method.    use butterfly or dgemm
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Nils Wedi + Mats Hamrud + George Modzynski
!
!     Modifications.
!     --------------
!        J.Hague : Oct 2012 DR_HOOK round calls to DGEMM:
!      F. Vana  05-Mar-2015  Support for single precision
!     ------------------------------------------------------------------
  
USE TPM_GEN         ,ONLY : LSYNC_TRANS
USE PARKIND_ECTRANS  ,ONLY : JPIM     ,JPRB,  JPRBT
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE TPM_DIM         ,ONLY : R, R_NDGNH,R_NSMAX
USE TPM_GEOMETRY    ,ONLY : G, G_NDGLU
USE TPM_FIELDS      ,ONLY : F, ZIA, &
     & ZAA,ZAS,LDZAA,LDZAS,TDZAA,TDZAS,&
     & IZBS,ILDZBA,ILDZBS,ITDZBA,ITDZBS,&
     & IZCS,IZCST,ILDZCA,ILDZCS,ITDZCA,ITDZCS,&
     & TDZAS, IF_FS_INV, ZAMAX, ZSMAX 
USE TPM_DISTR       ,ONLY : D,D_NUMP,D_MYMS
USE TPM_GEN, ONLY: NOUT
USE TPM_FLT
USE BUTTERFLY_ALG_MOD
USE CUDA_GEMM_BATCHED_MOD
USE, INTRINSIC :: ISO_C_BINDING
USE MPL_MODULE      ,ONLY : MPL_BARRIER
USE IEEE_ARITHMETIC

IMPLICIT NONE


!     DUMMY ARGUMENTS
INTEGER(KIND=JPIM)  :: KM
INTEGER(KIND=JPIM)  :: KMLOC
INTEGER(KIND=JPIM), INTENT(IN)  :: KFC
INTEGER(KIND=JPIM)  :: KIFC
INTEGER(KIND=JPIM)  :: KDGLU
INTEGER(KIND=JPIM), INTENT(IN)  :: KF_OUT_LT
REAL(KIND=JPRB),    INTENT(IN)  :: PIA(:,:,:)
REAL(KIND=JPRBT),    INTENT(OUT) :: PSOA1(:,:,:)
REAL(KIND=JPRBT),    INTENT(OUT) :: PAOA1(:,:,:)

!     LOCAL
INTEGER(KIND=JPIM) :: IA, ILA, ILS, IS, ISKIP, ISL, J1, IF, JGL,JK, J,JI, IRET
INTEGER(KIND=JPIM) :: ITHRESHOLD
REAL(KIND=JPRBT), ALLOCATABLE :: ZZBS(:,:,:), ZZCSTS(:,:,:)
REAL(KIND=JPRBT), ALLOCATABLE :: ZZBA(:,:,:), ZZCSTA(:,:,:)

INTEGER(KIND=JPIM) :: ISTAT

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  
!*       1.1      PREPARATIONS.
IF (LHOOK) CALL DR_HOOK('LE_DGEMM',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------

!*       1.1      PREPARATIONS.
IF (LSYNC_TRANS) THEN
  CALL MPL_BARRIER(CDSTRING='LEINV BARRIER')
ENDIF
CALL GSTATS(453,0)

ALLOCATE(ZZBA(ITDZBA,ILDZBA,D_NUMP))
ALLOCATE(ZZCSTA(ITDZCA,ILDZCA,D_NUMP))
ALLOCATE(ZZBS(ITDZBS,ILDZBS,D_NUMP))
ALLOCATE(ZZCSTS(ITDZCS,ILDZCS,D_NUMP))
!$ACC DATA CREATE(ZZBA,ZZCSTA,ZZBS,ZZCSTS)

!$ACC DATA COPYIN (S,S%ITHRESHOLD,S%LUSEFLT) &
!$ACC&      COPYIN (D,D_MYMS,R,G,G_NDGLU,D_NUMP,R_NDGNH,R_NSMAX) &
!$ACC&      PRESENT (ZAA,ZAS) &
!$ACC&      PRESENT (ZZBA,ZZBS,ZZCSTA,ZZCSTS) &
!$ACC&      PRESENT (PIA) &
!$ACC&      PRESENT (PSOA1,PAOA1,IZBS)

!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM) DEFAULT(NONE)
DO KMLOC=1,D_NUMP
   DO JGL=1,R_NDGNH
      DO J1=2,KFC,2
         
         KM = D_MYMS(KMLOC)
         IF(KM == 0)THEN
            PSOA1(J1,JGL,KMLOC) = 0.0_JPRBT
            PAOA1(J1,JGL,KMLOC) = 0.0_JPRBT
         END IF
      ENDDO
   ENDDO
   !end loop over wavenumber
END DO

!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,ISKIP,ILA,ILS,IA) DEFAULT(NONE)
DO KMLOC=1,D_NUMP
  DO J=1,(R_NSMAX+2)/2
    DO JK=1,KFC
   
      KM =  D_MYMS(KMLOC)
      IF (KM == 0) THEN
         ISKIP = 2
      ELSE
         ISKIP = 1
      ENDIF
      
      IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
         ILA = (R_NSMAX-KM+2)/2
         IF (J .LE. ILA) THEN
            IA  = 1+MOD(R_NSMAX-KM+2,2)
!!            IZBA((JK-1)/ISKIP+1,J,KMLOC)=PIA(JK,IA+1+(J-1)*2,KMLOC)*RRPELTMINV/ZAMAX((JK-1)/ISKIP+1,KMLOC)
            IZBS((JK-1)/ISKIP+1+(J-1+(KMLOC-1)*TDZAA)*IF_FS_INV)=PIA(JK,IA+1+(J-1)*2,KMLOC)
         ENDIF
      ENDIF
    ENDDO
  ENDDO
ENDDO

!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,ISKIP,ILA,ILS,IA) DEFAULT(NONE)
DO KMLOC=1,D_NUMP
  DO J=1,(R_NSMAX+2)/2
    DO JK=1,KFC
      KM =  D_MYMS(KMLOC)
      IF(KM == 0)THEN
         ISKIP = 2
      ELSE
         ISKIP = 1
      ENDIF
      
      IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
         ILA = (R_NSMAX-KM+2)/2
         IF (J .LE. ILA) THEN
            IA  = 1+MOD(R_NSMAX-KM+2,2)
            ZZBA((JK-1)/ISKIP+1,J,KMLOC)=PIA(JK,IA+1+(J-1)*2,KMLOC)
         ENDIF
      ENDIF
    ENDDO
  ENDDO
ENDDO

ITHRESHOLD=S%ITHRESHOLD

! operate on full arrays, where non-relavent entries have been set to zero
! call CUDA_DGEMM_BATCHED('N','N',LDZAA,TDZBA,TDZAA,1.0_JPRB,ZAA,LDZAA,TDZAA,ZBA,LDZBA,TDZBA,0._JPRB,ZCA,LDZCA,TDZCA,D_NUMP)
! Get C in transpose format to get better memory access patterns later
!C=A*B =>
! C^T=B^T*A^T


! OVERLOADED FOR SINGLE AND DOUBLE PRECISION
!$ACC HOST_DATA USE_DEVICE(ZAA,ZZBA,ZZCSTA)
CALL CUDA_GEMM_BATCHED( &
  & 'N', 'T', &
  & ITDZCA, ILDZCA, ILDZBA, &
  & 1.0_JPRBT, &
  & ZZBA, ITDZBA, ILDZBA,&
  & ZAA, LDZAA, TDZAA, &
  & 0._JPRBT, &
  & ZZCSTA, ITDZCA, ILDZCA, &
  & D_NUMP)
!$ACC END HOST_DATA

!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,KDGLU,ISKIP,ISL) DEFAULT(NONE)
DO KMLOC=1,D_NUMP
   DO JI=1,R_NDGNH
      DO JK=1,KFC
         KM = D_MYMS(KMLOC)
         KDGLU = MIN(R_NDGNH,G_NDGLU(KM))
         IF (JI .LE. KDGLU) then
            IF(KM == 0)THEN
               ISKIP = 2
            ELSE
               ISKIP = 1
            END IF
            
            ISL = MAX(R_NDGNH-G_NDGLU(KM)+1,1)
            IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
               PAOA1(JK,ISL+JI-1,KMLOC) = ZZCSTA((JK-1)/ISKIP+1,JI,KMLOC)
            END IF
         END IF
      ENDDO
   ENDDO
END DO

! 2. +++++++++++++ symmetric

!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,ISKIP,ILS,IS) DEFAULT(NONE)
DO KMLOC=1,D_NUMP
  DO J=1,(R_NSMAX+3)/2
    DO JK=1,KFC
      KM =  D_MYMS(KMLOC)
      IF(KM == 0)THEN
         ISKIP = 2
      ELSE
         ISKIP = 1
      ENDIF
      
      IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
         ILS = (R_NSMAX-KM+3)/2
         IF (J .LE. ILS) THEN
            IS  = 1+MOD(R_NSMAX-KM+1,2)
           ZZBS((JK-1)/ISKIP+1,J,KMLOC)=PIA(JK,IS+1+(J-1)*2,KMLOC)
         END IF
      END IF
    ENDDO
  ENDDO
ENDDO

!C=A*B =>
! C^T=B^T*A^T

!$ACC HOST_DATA USE_DEVICE(ZAS,ZZBS,ZZCSTS)
CALL CUDA_GEMM_BATCHED( &
  & 'N', 'T', &
  & ITDZCS, ILDZCS, ILDZBS, &
  & 1.0_JPRBT, &
  & ZZBS, ITDZBS, ILDZBS, &
  & ZAS, LDZAS, TDZAS, &
  & 0._JPRBT, &
  & ZZCSTS, ITDZCS, ILDZCS, &
  & D_NUMP)
!$ACC END HOST_DATA


!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,KDGLU,ISKIP,ISL) DEFAULT(NONE)
DO KMLOC=1,D_NUMP
   DO JI=1,R_NDGNH
      DO JK=1,KFC
         KM = D_MYMS(KMLOC)
         KDGLU = MIN(R_NDGNH,G_NDGLU(KM))
         IF (JI .LE. KDGLU) then
            IF(KM == 0)THEN
               ISKIP = 2
            ELSE
               ISKIP = 1
            END IF
            
            ISL = MAX(R_NDGNH-G_NDGLU(KM)+1,1)
            IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
               !PSOA1(JK,ISL+JI-1,KMLOC) = IZCST((JK-1)/ISKIP+1+(JI-1+(KMLOC-1)*R_NDGNH)*IF_FS_INV)
               PSOA1(JK,ISL+JI-1,KMLOC) = ZZCSTS((JK-1)/ISKIP+1,JI,KMLOC)
            END IF
         END IF
      ENDDO
   ENDDO
END DO

!$ACC END DATA

!$ACC END DATA
DEALLOCATE(ZZBS)
DEALLOCATE(ZZBA)
DEALLOCATE(ZZCSTS)
DEALLOCATE(ZZCSTA)

IF (LSYNC_TRANS) THEN
  CALL MPL_BARRIER(CDSTRING='LEINV BARRIER')
ENDIF
CALL GSTATS(453,1)

IF (LHOOK) CALL DR_HOOK('LE_DGEMM',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE LEINV
END MODULE LEINV_MOD
