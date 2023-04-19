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
      SUBROUTINE LEINV(KFC,KF_OUT_LT,PIA)

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
  
USE PARKIND_ECTRANS ,ONLY : JPIM ,JPRB, JPIB, JPRBT
USE YOMHOOK         ,ONLY : LHOOK,   DR_HOOK, JPHOOK
  
USE TPM_DIM         ,ONLY : R, R_NDGNH,R_NSMAX
USE TPM_GEOMETRY    ,ONLY : G, G_NDGLU
!!USE TPM_FIELDS    ,ONLY : F, ZIA, &
USE TPM_FIELDS      ,ONLY : ZSOA1,ZAOA1, &
 &                          ZAA,ZAS,LDZAA,LDZAS,TDZAA,TDZAS,&
 &                          IZBS,ILDZBA,ILDZBS,ITDZBA,ITDZBS,&
 &                          IZCS,IZCST,ILDZCA,ILDZCS,ITDZCA,ITDZCS,&
 &                          TDZAS,IF_FS_INV
  
USE TPM_DISTR       ,ONLY : D,D_NUMP,D_MYMS
      !
USE TPM_FLT
USE TPM_GEN         ,ONLY : NOUT ! Fpr nout
  
USE CUDA_GEMM_BATCHED_MOD
USE CUDA_DEVICE_MOD
  
USE OPENACC
USE ISO_C_BINDING
USE IEEE_ARITHMETIC

IMPLICIT NONE

      INTERFACE
         SUBROUTINE cudaProfilerStart() BIND(C,name='cudaProfilerStart')
           USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
           IMPLICIT NONE
         END SUBROUTINE cudaProfilerStart
      END INTERFACE

      INTERFACE
         SUBROUTINE cudaProfilerStop() BIND(C,name='cudaProfilerStop')
           USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT
           IMPLICIT NONE
         END SUBROUTINE cudaProfilerStop
      END INTERFACE
  
  
INTEGER(KIND=JPIM)  :: KM
INTEGER(KIND=JPIM)  :: KMLOC
INTEGER(KIND=JPIM), INTENT(IN)  :: KFC
INTEGER(KIND=JPIM)  :: KIFC
INTEGER(KIND=JPIM)  :: KDGLU
INTEGER(KIND=JPIM), INTENT(IN)  :: KF_OUT_LT
REAL(KIND=JPRB),    INTENT(IN)  :: PIA(:,:,:)
!REAL(KIND=JPRB),    INTENT(OUT) :: PSOA1(:,:,:)
!REAL(KIND=JPRB),    INTENT(OUT) :: PAOA1(:,:,:)

!     LOCAL
INTEGER(KIND=JPIM) :: IA, ILA, ILS, IS, ISKIP, ISL, J1, IF, JGL,JK, J,JI, IRET
INTEGER(KIND=JPIM) :: ITHRESHOLD

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  
!*       1.1      PREPARATIONS.
IF (LHOOK) CALL DR_HOOK('LE_DGEMM',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------
  

!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------

!*       1.1      PREPARATIONS.


#ifdef ACCGPU
!$ACC DATA COPYIN (S,S%ITHRESHOLD,S%LUSEFLT) &
!$ACC&     COPYIN (D,D_MYMS,R,G,G_NDGLU,D_NUMP,R_NDGNH,R_NSMAX) &
!$ACC&     PRESENT (ZAA,ZAS) &
!$ACC&     PRESENT (IZBS,IZCST) &
!$ACC&     PRESENT (ZSOA1,ZAOA1) &
!$ACC&     PRESENT (PIA)
#endif
#ifdef OMPGPU
!$OMP TARGET DATA MAP (TO:S,S%ITHRESHOLD,S%LUSEFLT) &
!$OMP&     MAP (TO:D,D_MYMS,R,G,G_NDGLU,D_NUMP,R_NDGNH,R_NSMAX) &
!$OMP&     MAP (ALLOC:ZAA,ZAS) &
!$OMP&     MAP (ALLOC:IZBS,IZCST) &
!$OMP&     MAP (ALLOC:ZSOA1,ZAOA1) &
!$OMP&     MAP (ALLOC:PIA)
#endif

#ifdef OMPGPU
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) PRIVATE(KM)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM)
#endif
DO KMLOC=1,D_NUMP
   DO JGL=1,R_NDGNH
      DO J1=2,KFC,2
         KM = D_MYMS(KMLOC)
         IF(KM == 0)THEN
                  ZSOA1(J1,JGL,KMLOC) = 0.0_JPRBT
                  ZAOA1(J1,JGL,KMLOC) = 0.0_JPRBT
         END IF
      ENDDO
   ENDDO
END DO

      ! 1. +++++++++++++ anti-symmetric
      
#ifdef OMPGPU
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) PRIVATE(KM,ISKIP,ILA,IA)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,ISKIP,ILA,IA)
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
                  IZBS((JK-1)/ISKIP+1+(J-1)*ITDZBA+(KMLOC-1)*ILDZBA*ITDZBA)=PIA(JK,IA+1+(J-1)*2,KMLOC)
         ENDIF
      ENDIF
    ENDDO
  ENDDO
ENDDO

ITHRESHOLD=S%ITHRESHOLD

! operate on full arrays, where non-relavent entries have been set to zero
! CALL HIP_DGEMM_BATCHED('N','N',LDZAA,TDZBA,TDZAA,1.0_JPRB,ZAA,LDZAA,TDZAA,ZBA,LDZBA,TDZBA,0._JPRB,ZCA,LDZCA,TDZCA,D_NUMP)
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
CALL HIP_GEMM_BATCHED('N','T',ITDZCA,ILDZCA,ILDZBA,1.0_JPRB,IZBS,ITDZBA,ILDZBA,&
 &                     ZAA,LDZAA,TDZAA,0._JPRB,IZCST,ITDZCA,ILDZCA,D_NUMP)
#ifdef OMPGPU
!$OMP END TARGET DATA
#endif
#ifdef ACCGPU
!$ACC END HOST_DATA
#endif

#ifdef OMPGPU
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) PRIVATE(KM,KDGLU,ISKIP,ISL)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,KDGLU,ISKIP,ISL)
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
               ZAOA1(JK,ISL+JI-1,KMLOC) = IZCST((JK-1)/ISKIP+1+(JI-1)*ITDZCA+(KMLOC-1)*ILDZCA*ITDZCA)
            END IF
         END IF
      ENDDO
   ENDDO
END DO

! 2. +++++++++++++ symmetric

#ifdef OMPGPU
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) PRIVATE(KM,ISKIP,ILS,IS)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,ISKIP,ILS,IS)
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
           IZBS((JK-1)/ISKIP+1+(J-1)*ITDZBS+(KMLOC-1)*ILDZBS*ITDZBS)=PIA(JK,IS+1+(J-1)*2,KMLOC)
         END IF
      END IF
    ENDDO
  ENDDO
ENDDO

!C=A*B =>
! C^T=B^T*A^T

#ifdef ACCGPU
!$ACC HOST_DATA USE_DEVICE(ZAS,IZBS,IZCST)
#endif
#ifdef OMPGPU
!$OMP TARGET DATA USE_DEVICE_PTR(ZAS,IZBS,IZCST)
#endif
CALL HIP_GEMM_BATCHED('N','T',ITDZCS,ILDZCS,ILDZBS,1.0_JPRB,IZBS,ITDZBS,ILDZBS,&
 &                     ZAS,LDZAS,TDZAS,0._JPRB,IZCST,ITDZCS,ILDZCS,D_NUMP)
#ifdef OMPGPU
!$OMP END TARGET DATA
#endif
#ifdef ACCGPU
!$ACC END HOST_DATA
#endif

!!istat = cuda_Synchronize()
#ifdef OMPGPU
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) PRIVATE(KM,KDGLU,ISKIP,ISL)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,KDGLU,ISKIP,ISL)
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
               ZSOA1(JK,ISL+JI-1,KMLOC) = IZCST((JK-1)/ISKIP+1+(JI-1)*ITDZCS+(KMLOC-1)*ITDZCS*ILDZCS)
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

IF (LHOOK) CALL DR_HOOK('LE_DGEMM',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE LEINV
END MODULE LEINV_MOD
