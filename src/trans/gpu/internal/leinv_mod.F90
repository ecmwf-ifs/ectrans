! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LEINV_MOD
CONTAINS
SUBROUTINE LEINV(KFC,KSTA,KF_OUT_LT,PAOA1,PSOA1)

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
  
USE PARKIND_ECTRANS ,ONLY : JPIM     ,JPRB,  JPRBT
USE YOMHOOK         ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE TPM_DIM         ,ONLY : R, R_NDGNH,R_NSMAX
USE TPM_GEOMETRY    ,ONLY : G, G_NDGLU
USE TPM_FIELDS      ,ONLY : F, ZIA, &
 &                          ZAA,ZAS,LDZAA,LDZAS,TDZAA,TDZAS,&
 &                          IZBS,ILDZBA,ILDZBS,ITDZBA,ITDZBS,&
 &                          IZCS,IZCST,ILDZCA,ILDZCS,ITDZCA,ITDZCS,&
 &                          TDZAS, IF_FS_INV, ZAMAX, ZSMAX 
USE TPM_DISTR       ,ONLY : D_NUMP,D_MYMS, MYPROC
USE TPM_GEN         ,ONLY : NOUT
USE TPM_FLT
USE CUBLAS_MOD       ,ONLY : CUDA_SGEMM_BATCHED
! issue ? address error
!USE CUDA_GEMM_BATCHED_MOD
USE, INTRINSIC :: ISO_C_BINDING
USE IEEE_ARITHMETIC

IMPLICIT NONE


!     DUMMY ARGUMENTS
INTEGER(KIND=JPIM)  :: KM
INTEGER(KIND=JPIM)  :: KMLOC
INTEGER(KIND=JPIM), INTENT(IN)  :: KSTA
INTEGER(KIND=JPIM), INTENT(IN)  :: KFC
INTEGER(KIND=JPIM)  :: KIFC
INTEGER(KIND=JPIM)  :: KDGLU
INTEGER(KIND=JPIM), INTENT(IN)  :: KF_OUT_LT
!REAL(KIND=JPRB),    INTENT(IN)  :: PIA(:,:,:)
REAL(KIND=JPRBT),    INTENT(OUT) :: PSOA1(:,:,:)
REAL(KIND=JPRBT),    INTENT(OUT) :: PAOA1(:,:,:)

!     LOCAL
INTEGER(KIND=JPIM) :: IA, ILA, ILS, IS, ISKIP, ISL, J1, IF, JGL,JK, J,JI, IRET
INTEGER(KIND=JPIM) :: ITHRESHOLD

INTEGER(KIND=JPIM) :: ISTAT

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
!*       1.1      PREPARATIONS.
IF (LHOOK) CALL DR_HOOK('LE_DGEMM',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------

!*       1.1      PREPARATIONS.
#ifdef ACCGPU
!$ACC DATA COPYIN (S,S%ITHRESHOLD,S%LUSEFLT) &
!$ACC&      COPYIN (D_MYMS,G_NDGLU,D_NUMP,R_NDGNH,R_NSMAX, KFC, KSTA) &
!$ACC&      COPYIN (IF_FS_INV) &
!$ACC&      PRESENT (ZIA,ZAA,ZAS) &
!$ACC&      PRESENT (IZCST,IZBS) &
!!$ACC&      COPYIN (PIA) &
!$ACC&      COPYOUT (PSOA1,PAOA1)
#endif
#ifdef OMPGPU
!$OMP TARGET DATA MAP(TO:S,S%ITHRESHOLD,S%LUSEFLT,D_MYMS,G_NDGLU,D_NUMP,R_NDGNH,R_NSMAX) &
!$OMP&      MAP(ALLOC:ZAA,ZAS,IZCST,ZIA,PSOA1,PAOA1,IZBS)
#endif

#ifdef OMPGPU
!$OMP TARGET PARALLEL DO COLLAPSE(3) PRIVATE(KM) &
!$OMP&      SHARED(D_MYMS,PSOA1,PAOA1,KFC,R_NDGNH,D_NUMP) DEFAULT(NONE)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,KMLOC,JGL) DEFAULT(NONE) &
!$ACC&      PRESENT(D_MYMS,PSOA1,PAOA1,KFC,KSTA,R_NDGNH,D_NUMP)
#endif
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

#ifdef OMPGPU
!$OMP TARGET PARALLEL DO COLLAPSE(3) PRIVATE(KM,ISKIP,ILA,ILS,IA) &
!$OMP&      SHARED(D_NUMP,R_NSMAX,KFC,D_MYMS,IZBS,TDZAA,IF_FS_INV,ZIA) DEFAULT(NONE)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KMLOC,J,JK,KM,ISKIP,ILA,ILS,IA) DEFAULT(NONE) &
!$ACC&      PRESENT(D_NUMP,R_NSMAX,KFC,KSTA,D_MYMS,IZBS,IF_FS_INV,ZIA,TDZAA)
#endif
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
            IZBS((JK-1)/ISKIP+1+(J-1+(KMLOC-1)*TDZAA)*IF_FS_INV)=ZIA(KSTA+JK-1,IA+1+(J-1)*2,KMLOC)
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
#ifdef ACCGPU
!$ACC HOST_DATA USE_DEVICE(ZAA,IZBS,IZCST)
#endif
#ifdef OMPGPU
!$OMP TARGET DATA USE_DEVICE_PTR(ZAA,IZBS,IZCST)
#endif
CALL CUDA_SGEMM_BATCHED( &
  & 'N', 'T', &
  & ITDZCA, ILDZCA, ILDZBA, &
  & 1.0_JPRBT, &
  & IZBS, ITDZBA, ILDZBA,&
  & ZAA, LDZAA, TDZAA, &
  & 0._JPRBT, &
  & IZCST, ITDZCA, ILDZCA, &
  & D_NUMP)
#ifdef OMPGPU
!$OMP END TARGET DATA
#endif
#ifdef ACCGPU
!$ACC END HOST_DATA
#endif

#ifdef OMPGPU
!$OMP TARGET PARALLEL DO COLLAPSE(3) PRIVATE(KM,KDGLU,ISKIP,ISL) &
!$OMP&      SHARED(D_NUMP,R_NDGNH,KFC,D_MYMS,G_NDGLU,PAOA1,IZCST,ITDZCA,ILDZCA) DEFAULT(NONE)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,KDGLU,ISKIP,ISL) DEFAULT(NONE) &
!$ACC&      PRESENT(D_NUMP,R_NDGNH,KFC,D_MYMS,G_NDGLU,PAOA1,IZCST,ITDZCA,ILDZCA)
#endif
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
               PAOA1(JK,ISL+JI-1,KMLOC) = IZCST((JK-1)/ISKIP+1+(JI-1)*ITDZCA+(KMLOC-1)*ILDZCA*ITDZCA)
            END IF
         END IF
      ENDDO
   ENDDO
END DO

! 2. +++++++++++++ symmetric

#ifdef OMPGPU
!$OMP TARGET PARALLEL DO COLLAPSE(3) PRIVATE(KM,ISKIP,ILS,IS) &
!$OMP&      SHARED(D_NUMP,R_NSMAX,KFC,D_MYMS,IZBS,ITDZBS,ILDZBS,ZIA) DEFAULT(NONE)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,ISKIP,ILS,IS) DEFAULT(NONE) &
!$ACC&      PRESENT(D_NUMP,R_NSMAX,KFC,KSTA,D_MYMS,IZBS,ITDZBS,ILDZBS,ZIA)
#endif
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
           IZBS((JK-1)/ISKIP+1+(J-1)*ITDZBS+(KMLOC-1)*ILDZBS*ITDZBS)=ZIA(KSTA+JK-1,IS+1+(J-1)*2,KMLOC)
         END IF
      END IF
    ENDDO
  ENDDO
ENDDO

!C=A*B =>
! C^T=B^T*A^T

#ifdef OMPGPU
!$OMP TARGET DATA USE_DEVICE_PTR(ZAS,IZBS,IZCST)
#endif
#ifdef ACCGPU
!$ACC HOST_DATA USE_DEVICE(ZAS,IZBS,IZCST)
#endif
CALL CUDA_SGEMM_BATCHED( &
 & 'N', 'T', &
 & ITDZCS, ILDZCS, ILDZBS, &
 & 1.0_JPRBT, &
 & IZBS, ITDZBS, ILDZBS, &
 & ZAS, LDZAS, TDZAS, &
 & 0._JPRBT, &
 & IZCST, ITDZCS, ILDZCS, &
 & D_NUMP)
#ifdef OMPGPU
!$OMP END TARGET DATA
#endif
#ifdef ACCGPU
!$ACC END HOST_DATA
#endif


#ifdef OMPGPU
!$OMP TARGET PARALLEL DO COLLAPSE(3) PRIVATE(KM,KDGLU,ISKIP,ISL) &
!$OMP&      SHARED(D_NUMP,R_NDGNH,KFC,D_MYMS,G_NDGLU,PSOA1,IZCST,ITDZCS,ILDZCS) DEFAULT(NONE)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,KDGLU,ISKIP,ISL) DEFAULT(NONE) &
!$ACC&      PRESENT(D_NUMP,R_NDGNH,KFC,D_MYMS,G_NDGLU,PSOA1,IZCST,ITDZCS,ILDZCS)
#endif
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
               PSOA1(JK,ISL+JI-1,KMLOC) = IZCST((JK-1)/ISKIP+1+(JI-1)*ITDZCS+(KMLOC-1)*ITDZCS*ILDZCS)
            END IF
         END IF
      ENDDO
   ENDDO
END DO

#ifdef OMPGPU
!$OMP END TARGET DATA
#endif
#ifdef ACCGPU
!$ACC END DATA
#endif

!*       1.       PERFORM LEGENDRE TRANFORM.

IF (LHOOK) CALL DR_HOOK('LE_DGEMM',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE LEINV
END MODULE LEINV_MOD
