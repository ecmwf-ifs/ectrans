! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LEDIR_MOD
CONTAINS
SUBROUTINE LEDIR(KF_FS,KLED2,PAIA,POA1,KMODE)

!**** *LEDIR* - Direct Legendre transform.

!     Purpose.
!     --------
!        Direct Legendre tranform of state variables.

!**   Interface.
!     ----------
!        CALL LEDIR(...)

!        Explicit arguments :  KM - zonal wavenumber
!        --------------------  KFC - number of field to transform
!                              PAIA - antisymmetric part of Fourier
!                              fields for zonal wavenumber KM
!                              PSIA - symmetric part of Fourier
!                              fields for zonal wavenumber KM
!                              POA1 -  spectral
!                              fields for zonal wavenumber KM

!        Implicit arguments :  None.
!        --------------------

!     Method.
!     -------   use butterfly or dgemm

!     Externals.   
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!          Nils Wedi + Mats Hamrud + George Modzynski

!     Modifications.
!     --------------
!        J.Hague : Oct 2012 DR_HOOK round calls to DGEMM:
!      F. Vana  05-Mar-2015  Support for single precision
!     ------------------------------------------------------------------

USE PARKIND_ECTRANS  ,ONLY : JPIM ,JPIB    ,JPRB,  JPRBT
USE YOMHOOK          ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE TPM_DIM          ,ONLY : R_NDGNH,R_NSMAX,R_NTMAX
USE TPM_GEOMETRY     ,ONLY : G_NDGLU
USE TPM_FIELDS       ,ONLY : F_RW, &
 &                           ZAA,ZAS,LDZAA,LDZAS,TDZAA,TDZAS,&
 &                           DZBST,DLDZBA,DLDZBS,DTDZBA,DTDZBS,&
 &                           DZCST,DZCAT,DLDZCA,DLDZCS,DTDZCA,DTDZCS,&
 &                           ZAMAX, ZSMAX,&
 &                           IF_FS_DIR,ZAA0,DZBST0,DZCAT0,ZAS0,DZCST0,KMLOC0
USE TPM_DISTR
USE TPM_GEN          ,ONLY : NOUT
USE TPM_FLT
USE BUTTERFLY_ALG_MOD
USE HIPBLAS_MOD       ,ONLY : HIP_DGEMM_BATCHED, HIP_SGEMM_BATCHED
!USE HIP_GEMM_BATCHED_MOD!!, ONLY: HIP_TCGEMM_BATCHED, HIP_GEMM_BATCHED
USE, INTRINSIC :: ISO_C_BINDING
USE IEEE_ARITHMETIC

IMPLICIT NONE


!     DUMMY ARGUMENTS
INTEGER(KIND=JPIM)  :: KM
INTEGER(KIND=JPIM)  :: KMLOC
INTEGER(KIND=JPIM)  :: KFC
INTEGER(KIND=JPIM)  :: KIFC
INTEGER(KIND=JPIM) :: KDGLU
INTEGER(KIND=JPIM), INTENT(IN)  :: KF_FS
INTEGER(KIND=JPIM), INTENT(IN)  :: KLED2
INTEGER(KIND=JPIM), INTENT(IN)  :: KMODE

REAL(KIND=JPRBT),    INTENT(IN)  :: PAIA(:,:,:)
!REAL(KIND=JPRBT),    INTENT(IN)  :: PSIA(:,:,:),   PAIA(:,:,:)
REAL(KIND=JPRBT),    INTENT(OUT) :: POA1(:,:,:)

!     LOCAL VARIABLES
INTEGER(KIND=JPIM) :: IA, ILA, ILS, IS, ISKIP, ISL, IF, J, JK, IRET
INTEGER(KIND=JPIM) :: ITHRESHOLD
REAL(KIND=JPRB) :: RRPELTMDIR = 100.0_JPRB
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

INTEGER :: ISTAT

IF (LHOOK) CALL DR_HOOK('LE_DGEMM',0,ZHOOK_HANDLE)

KFC = 2*KF_FS
KIFC = KFC

#ifdef ACCGPU
!$ACC DATA &
!$ACC& PRESENT(F_RW) &
!$ACC& PRESENT(D_NUMP,D_MYMS,R_NDGNH,G_NDGLU,R_NSMAX,R_NTMAX) &
!$ACC& PRESENT(PAIA) &
!$ACC& PRESENT(ZAA,ZAS,DZBST,DZCST,DZCAT) &
!$ACC& PRESENT(POA1,dzbst0,dzcat0,dzbst0,dzcst0)
#endif
#ifdef OMPGPU
!$OMP TARGET DATA &
!$OMP& MAP(TO:F_RW) &
!$OMP& MAP(TO:D_NUMP,D_MYMS,R_NDGNH,G_NDGLU,R_NSMAX,R_NTMAX) &
!$OMP& MAP(PRESENT,ALLOC:PAIA) &
!$OMP& MAP(PRESENT,ALLOC:ZAA,ZAS,DZBST,DZCST,DZCAT) &
!$OMP& MAP(PRESENT,ALLOC:POA1,dzbst0,dzcat0,dzbst0,dzcst0)
#endif


!! Initialize rescaling arrays to zero
!!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2)
!!$ACC PARALLEL LOOP COLLAPSE(2)
!DO KMLOC=1,SIZE(ZAMAX,2)
!  DO JK=1,SIZE(ZAMAX,1)
!    ZAMAX(JK,KMLOC) = 0.0_JPRBT
!    ZSMAX(JK,KMLOC) = 0.0_JPRBT
!  ENDDO
!ENDDO


! anti-symmetric

IF ( KMODE == -1 ) THEN

#ifdef OMPGPU
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) PRIVATE(KM,KDGLU,ISL,ISKIP) DEFAULT(NONE) &
!$OMP&     SHARED(D_NUMP,R_NDGNH,KFC,D_MYMS,G_NDGLU,DZBST,DLDZBA,DTDZBA,PAIA,F_RW)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,KDGLU,ISL,ISKIP) DEFAULT(NONE) &
!$ACC&     COPYIN(KFC,DZBST,DLDZBA,DTDZBA) &
!$ACC&     PRESENT(D_NUMP,R_NDGNH,D_MYMS,G_NDGLU,PAIA,F_RW)
#endif
DO KMLOC=1,D_NUMP
   DO J=1,R_NDGNH
      DO JK=1,KFC
         
         KM = D_MYMS(KMLOC)   
         KDGLU = MIN(R_NDGNH,G_NDGLU(KM))
         IF (J .LE. KDGLU) THEN
            ISL = MAX(R_NDGNH-G_NDGLU(KM)+1,1)
            
            IF(KM == 0)THEN
               ISKIP = 2
            ELSE
               ISKIP = 1
            ENDIF
            IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
               DZBST((JK-1)/ISKIP+1+(J-1+(KMLOC-1)*DLDZBA)*DTDZBA)=PAIA(JK,ISL+J-1,KMLOC)*F_RW(ISL+J-1)
            END IF
         END IF
      ENDDO
   ENDDO
END DO


! Get C in transpose format to get better memory access patterns later
!C=A*B =>
! C^T=B^T*A^T
#ifdef OMPGPU
!$OMP TARGET DATA USE_DEVICE_PTR(ZAA,DZBST,DZCAT)
#endif
#ifdef ACCGPU
!$ACC HOST_DATA USE_DEVICE(ZAA,DZBST,DZCAT)
#endif
CALL HIP_SGEMM_BATCHED( &
  & 'N', 'N', &
  & DTDZBA, TDZAA, DLDZBA, &
  & 1.0_JPRBT, &
  & DZBST, DTDZBA, DLDZBA, &
  & ZAA, LDZAA, TDZAA, &
  & 0._JPRBT, &
  & DZCAT, DTDZCA, DLDZCA, &
  & D_NUMP)
#ifdef OMPGPU
!$OMP END TARGET DATA
#endif
#ifdef ACCGPU
!$ACC END HOST_DATA
#endif

#ifdef OMPGPU
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) PRIVATE(KM,ISKIP,ILA,IA,ILS) DEFAULT(NONE) &
!$OMP&    SHARED(D_NUMP,R_NTMAX,KFC,D_MYMS,POA1,DZCAT,DLDZCA,DTDZCA)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,ISKIP,ILA,IA,ILS) DEFAULT(NONE) &
!$ACC&    COPYIN(KFC,DLDZCA,DTDZCA) &
!$ACC&    PRESENT(D_NUMP,R_NTMAX,D_MYMS,POA1,DZCAT)
#endif
DO KMLOC=1,D_NUMP
   DO J=1,(R_NTMAX+2)/2
      DO JK=1,KFC

         KM = D_MYMS(KMLOC)
         IF(KM == 0)THEN
            ISKIP = 2
         ELSE
            ISKIP = 1
         ENDIF
         
         IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
            ILA = (R_NTMAX-KM+2)/2
            IA  = 1+MOD(R_NTMAX-KM+2,2)
            IF (J .LE. ILA) THEN
               POA1(JK,IA+(J-1)*2,KMLOC) = DZCAT((JK-1)/ISKIP+1+(J-1+(KMLOC-1)*DLDZCA)*DTDZCA)
            END IF
         END IF
      ENDDO
   ENDDO
ENDDO

! compute m=0 in double precision:
IF(KMLOC0 > 0) THEN
   PRINT*,'computing m=0 in double precision'
   ISKIP = 2

#ifdef OMPGPU
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(KDGLU,ISL) DEFAULT(NONE) &
  !$OMP&     SHARED(R_NDGNH,KFC,G_NDGLU,DZBST0,DTDZBA,PAIA,KMLOC0,F_RW,ISKIP)
#endif
#ifdef ACCGPU
  !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KDGLU,ISL) DEFAULT(NONE) &
  !$ACC&     COPYIN(KFC,DTDZBA,KMLOC0,ISKIP) &
  !$ACC&     PRESENT(R_NDGNH,G_NDGLU,DZBST0,PAIA,F_RW)
#endif
   DO J=1,R_NDGNH
    DO JK=1,KFC

         KDGLU = MIN(R_NDGNH,G_NDGLU(0))
         IF (J .LE. KDGLU) THEN
            ISL = MAX(R_NDGNH-G_NDGLU(0)+1,1)
            IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
               DZBST0((JK-1)/ISKIP+1+(J-1)*DTDZBA)=PAIA(JK,ISL+J-1,KMLOC0)*F_RW(ISL+J-1)
            END IF
         END IF
    ENDDO
  ENDDO


  ! Get C in transpose format to get better memory access patterns later
  !C=A*B =>
  ! C^T=B^T*A^T

#ifdef OMPGPU
  !$OMP TARGET DATA USE_DEVICE_PTR(ZAA0,DZBST0,DZCAT0)
#endif
#ifdef ACCGPU
  !$ACC HOST_DATA USE_DEVICE(ZAA0,DZBST0,DZCAT0)
#endif
  CALL HIP_DGEMM_BATCHED( &
       & 'N','N', &
       & DTDZBA, TDZAA, DLDZBA, &
       & 1.0_JPRD, &
       & DZBST0,DTDZBA, DLDZBA,&
       & ZAA0,LDZAA, TDZAA, &
       & 0._JPRD, &
       & DZCAT0,DTDZCA, DLDZCA, &
       & 1)
  !CALL HIP_DGEMM('N','N',DTDZBA,TDZAA,DLDZBA,1.0_JPRD,DZBST0,DTDZBA,&
  !      &ZAA0,LDZAA,0._JPRD,DZCAT0,DTDZCA)
#ifdef OMPGPU
  !$OMP END TARGET DATA
#endif
#ifdef ACCGPU
  !$ACC END HOST_DATA
#endif

#ifdef OMPGPU
   !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(ILA,IA,ILS) DEFAULT(NONE) &
   !$OMP&    SHARED(R_NTMAX,KFC,POA1,DZCAT0,DTDZCA,KMLOC0,ISKIP)
#endif
#ifdef ACCGPU
   !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(ILA,IA,ILS) DEFAULT(NONE) &
   !$ACC&    COPYIN(KFC,DTDZCA,KMLOC0,ISKIP) &
   !$ACC&    PRESENT(R_NTMAX,POA1,DZCAT0)
#endif
   DO J=1,(R_NTMAX+2)/2
   DO JK=1,KFC
         IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
            ILA = (R_NTMAX+2)/2
            IA  = 1+MOD(R_NTMAX+2,2)
            IF (J .LE. ILA) THEN
               POA1(JK,IA+(J-1)*2,KMLOC0) = DZCAT0((JK-1)/ISKIP+1+(J-1)*DTDZCA)
            END IF
         END IF
   ENDDO
ENDDO
ENDIF

ELSE

! symmetric

#ifdef OMPGPU
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) PRIVATE(KM,KDGLU,ISL,ISKIP) DEFAULT(NONE) &
!$OMP&     SHARED(D_NUMP,R_NDGNH,KFC,D_MYMS,G_NDGLU,DZBST,DLDZBS,DTDZBS,PAIA,F_RW)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,KDGLU,ISL,ISKIP) DEFAULT(NONE) &
!$ACC &    COPYIN(KFC,DLDZBS,DTDZBS) &
!$ACC &    PRESENT(D_NUMP,R_NDGNH,D_MYMS,G_NDGLU,DZBST,PAIA,F_RW)
#endif
DO KMLOC=1,D_NUMP
   DO J=1,R_NDGNH   
      DO JK=1,KFC
         KM = D_MYMS(KMLOC)   
         KDGLU = MIN(R_NDGNH,G_NDGLU(KM))
         IF (J .LE. KDGLU) THEN
            ISL = MAX(R_NDGNH-G_NDGLU(KM)+1,1)
            
            IF(KM == 0)THEN
               ISKIP = 2
            ELSE
               ISKIP = 1
            ENDIF
            IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
               DZBST((JK-1)/ISKIP+1+(J-1+(KMLOC-1)*DLDZBS)*DTDZBS)=PAIA(JK,ISL+J-1,KMLOC)*F_RW(ISL+J-1)
            END IF
         END IF
      ENDDO
   ENDDO
END DO

! Get C in transpose format to get better memory access patterns later
!C=A*B =>
! C^T=B^T*A^T
#ifdef OMPGPU
!$OMP TARGET DATA USE_DEVICE_PTR(ZAS,DZBST,DZCST)
#endif
#ifdef ACCGPU
!$ACC HOST_DATA USE_DEVICE(ZAS,DZBST,DZCST)
#endif
CALL HIP_SGEMM_BATCHED( &
  & 'N', 'N', &
  & DTDZBS, TDZAS, DLDZBS, &
  & 1.0_JPRBT, &
  & DZBST, DTDZBS, DLDZBS, &
  & ZAS, LDZAS, TDZAS, &
  & 0._JPRBT, &
  & DZCST, DTDZCS, DLDZCS, &
  & D_NUMP)
#ifdef OMPGPU
!$OMP END TARGET DATA
#endif
#ifdef ACCGPU
!$ACC END HOST_DATA
#endif

#ifdef OMPGPU
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) PRIVATE(KM,ISKIP,ILA,IS,ILS) DEFAULT(NONE) &
!$OMP&    SHARED(D_NUMP,R_NTMAX,KFC,D_MYMS,POA1,DZCST,DLDZCS,DTDZCS)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,ISKIP,ILA,ILS,IS) DEFAULT(NONE) &
!$ACC&    COPYIN(KFC,DLDZCS,DTDZCS) &
!$ACC&    PRESENT(D_NUMP,R_NTMAX,D_MYMS,POA1,DZCST)
#endif
  DO KMLOC=1,D_NUMP
   DO J=1,(R_NTMAX+3)/2
      DO JK=1,KFC

         KM = D_MYMS(KMLOC)
         IF(KM == 0)THEN
            ISKIP = 2
         ELSE
            ISKIP = 1
         ENDIF
         
         IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
            ILS = (R_NTMAX-KM+3)/2
            IF (J .LE. ILS) THEN
               IS  = 1+MOD(R_NTMAX-KM+1,2)
               POA1(JK,IS+(J-1)*2,KMLOC) = DZCST((JK-1)/ISKIP+1+(J-1+(KMLOC-1)*DLDZCS)*DTDZCS)            
            END IF
         END IF
      ENDDO
   ENDDO
ENDDO

IF(KMLOC0 > 0) THEN
   ISKIP=2
#ifdef OMPGPU
   !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(KDGLU,ISL) DEFAULT(NONE) &
   !$OMP&     SHARED(R_NDGNH,KFC,G_NDGLU,DZBST0,DTDZBS,PAIA,KMLOC0,F_RW,ISKIP)
#endif
#ifdef ACCGPU
   !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KDGLU,ISL) DEFAULT(NONE) &
   !$ACC&     COPYIN(KFC,DTDZBS,KMLOC0,ISKIP) &
   !$ACC&     PRESENT(R_NDGNH,G_NDGLU,DZBST0,PAIA,F_RW)
#endif
   DO J=1,R_NDGNH   
      DO JK=1,KFC
         KDGLU = MIN(R_NDGNH,G_NDGLU(0))
         IF (J .LE. KDGLU) THEN
            ISL = MAX(R_NDGNH-G_NDGLU(0)+1,1)
           IF (MOD((JK-1),ISKIP) .eq. 0) THEN
               DZBST0((JK-1)/ISKIP+1+(J-1)*DTDZBS)=PAIA(JK,ISL+J-1,KMLOC0)*F_RW(ISL+J-1)
           END IF
         END IF
      ENDDO
   ENDDO

      ! Get C in transpose format to get better memory access patterns later
      !C=A*B =>
      ! C^T=B^T*A^T

#ifdef ACCGPU
      !$ACC HOST_DATA USE_DEVICE(ZAS0,DZBST0,DZCST0)
#endif
#ifdef OMPGPU
      !$OMP TARGET DATA USE_DEVICE_PTR(ZAS0,DZBST0,DZCST0)
#endif
      CALL HIP_DGEMM_BATCHED('N','N',&
 &                            DTDZBS,TDZAS,DLDZBS,&
 &                            1.0_JPRD,DZBST0,DTDZBS,DLDZBS,&
 &                            ZAS0,LDZAS,TDZAS,&
 &                            0._JPRD,DZCST0,DTDZCS,DLDZCS,1)
#ifdef OMPGPU
      !$OMP END TARGET DATA
#endif
#ifdef ACCGPU
      !$ACC END HOST_DATA
#endif

#ifdef OMPGPU
   !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(ILA,ILS,IS) DEFAULT(NONE) &
   !$OMP&    SHARED(R_NTMAX,KFC,POA1,DZCST0,DTDZCS,KMLOC0,ISKIP)
#endif
#ifdef ACCGPU
   !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(ILA,ILS,IS) DEFAULT(NONE) &
   !$ACC&    COPYIN(KFC,DTDZCS,KMLOC0,ISKIP) &
   !$ACC&    PRESENT(R_NTMAX,POA1,DZCST0)
#endif
   DO J=1,(R_NTMAX+3)/2
      DO JK=1,KFC
         IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
            ILS = (R_NTMAX+3)/2
            IF (J .LE. ILS) THEN
               IS  = 1+MOD(R_NTMAX+1,2)
               POA1(JK,IS+(J-1)*2,KMLOC0) = DZCST0((JK-1)/ISKIP+1+(J-1)*DTDZCS)
            end if
         end if
      ENDDO
   ENDDO

ENDIF

ENDIF

#ifdef OMPGPU
!$OMP END TARGET DATA
#endif
#ifdef ACCGPU
!$ACC END DATA
#endif


IF (LHOOK) CALL DR_HOOK('LE_DGEMM',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE LEDIR
END MODULE LEDIR_MOD
