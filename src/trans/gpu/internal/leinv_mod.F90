! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
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
SUBROUTINE LEINV(KF_LT,PIA)

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
USE TPM_DIM         ,ONLY : R_NDGNH,R_NSMAX, R_NDGL
USE TPM_GEOMETRY    ,ONLY : G_NDGLU
USE TPM_FIELDS      ,ONLY : ZAA,ZAS,TDZAA,TDZAS
USE TPM_DISTR       ,ONLY : D_NUMP,D_MYMS, D_NPROCL,D_NPNTGTB1,D_NSTAGT1B
USE TPM_GEN         ,ONLY : NOUT, LSYNC_TRANS
USE TPM_TRANS       ,ONLY : FOUBUF_IN
USE TPM_FLT
USE HICBLAS_MOD     ,ONLY : HIP_SGEMM_BATCHED, HIP_DGEMM_BATCHED
USE MPL_MODULE      ,ONLY : MPL_BARRIER
#ifdef TRANS_SINGLE
#define HIP_GEMM_BATCHED HIP_SGEMM_BATCHED
#else
#define HIP_GEMM_BATCHED HIP_DGEMM_BATCHED
#endif
USE, INTRINSIC :: ISO_C_BINDING
USE IEEE_ARITHMETIC

IMPLICIT NONE


!     DUMMY ARGUMENTS
INTEGER(KIND=JPIM), INTENT(IN)  :: KF_LT
REAL(KIND=JPRB),    INTENT(IN)  :: PIA(:,:,:)

!     LOCAL
INTEGER(KIND=JPIM)  :: KM
INTEGER(KIND=JPIM)  :: KMLOC
INTEGER(KIND=JPIM)  :: KIFC
INTEGER(KIND=JPIM) :: IA, IS, ISL, J1, JGL, JK, J, IGLS, ISTAS, OFFSET1, OFFSET2
REAL(KIND=JPRBT), ALLOCATABLE :: ZINP(:), ZOUTS(:), ZOUTA(:)
REAL(KIND=JPRBT) :: ZAOA, ZSOA


INTEGER(KIND=JPIM) :: ISTAT

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('LE_DGEMM',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------


!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------

!*       1.1      PREPARATIONS.

IF (LSYNC_TRANS) THEN
  CALL MPL_BARRIER(CDSTRING='LEINV BARRIER')
ENDIF
CALL GSTATS(453,0)
ALLOCATE(ZINP(2*KF_LT*TDZAS*D_NUMP))
ALLOCATE(ZOUTS(2*KF_LT*R_NDGNH*D_NUMP))
ALLOCATE(ZOUTA(2*KF_LT*R_NDGNH*D_NUMP))
 
#ifdef ACCGPU
!$ACC DATA  COPYIN (D_MYMS,G_NDGLU,D_NUMP,R_NDGNH,R_NSMAX) &
!$ACC&      CREATE  (ZINP,ZOUTS,ZOUTA) &
!$ACC&      PRESENT (ZAA,ZAS,PIA,FOUBUF_IN)  &
!$ACC&      PRESENT(D_MYMS,D_NPROCL,D_NSTAGT1B,D_NPNTGTB1,G_NDGLU)
#endif
#ifdef OMPGPU
!$OMP TARGET DATA MAP(TO:D_MYMS,G_NDGLU,D_NUMP,R_NDGNH,R_NSMAX,ITDZBA, ITDZBS, ITDZCA, ITDZCS) &
!$OMP&      MAP(PRESENT,ALLOC:ZAA,ZAS,IZCST,ZIA,PSOA1,PAOA1,IZBS)
#endif

!IF KM=0 and NSMAX is 6:
!    IA=1
!    DO=1,6/2+1 ... 1..4
!       PIA_2=1+1+(J-1)*2 ...2+(0..3)*2 .... 2,4,6,8
!IF KM=0 and NSMAX is 7:
!    IA=2
!    DO=1,7/2+1 ... 1..4
!       PIA_2=2+1+(1..4-1)*2 ...3+(0..3)*2 .... 3,5,7,9

#ifdef OMPGPU
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) PRIVATE(KM,ISKIP,ILA,ILS,IA) &
!$OMP&      SHARED(D_NUMP,R_NSMAX,KSTA,D_MYMS,IZBS,TDZAA,IF_FS_INV,ZIA) DEFAULT(NONE)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM,IA,J) DEFAULT(NONE) &
!$ACC&    FIRSTPRIVATE(KF_LT,TDZAA)
#endif
DO KMLOC=1,D_NUMP
  DO JK=1,2*KF_LT
    KM =  D_MYMS(KMLOC)
    IA  = 1+MOD(R_NSMAX-KM+2,2)
    IF(KM /= 0)THEN
      !$ACC LOOP SEQ
      DO J=1,(R_NSMAX-KM+2)/2
        ZINP(JK+(J-1+(KMLOC-1)*TDZAA)*2*KF_LT)=PIA(JK,IA+1+(J-1)*2,KMLOC)
      ENDDO
    ELSEIF (MOD((JK+1),2) .EQ. 0) THEN
      !$ACC LOOP SEQ
      DO J=1,(R_NSMAX+2)/2
        ZINP((JK-1)/2+1+(J-1+(KMLOC-1)*TDZAA)*2*KF_LT)=PIA(JK,IA+1+(J-1)*2,KMLOC)
      ENDDO
    ENDIF
  ENDDO
ENDDO

! operate on full arrays, where non-relavent entries have been set to zero
! Get C in transpose format to get better memory access patterns later
!C=A*B =>
! C^T=B^T*A^T

! OVERLOADED FOR SINGLE AND DOUBLE PRECISION
#ifdef ACCGPU
!$ACC HOST_DATA USE_DEVICE(ZAA,ZINP,ZOUTA)
#endif
#ifdef OMPGPU
!$OMP TARGET DATA USE_DEVICE_PTR(ZAA,ZINP,ZOUTA)
#endif
CALL HIP_GEMM_BATCHED( &
  & 'N', 'T', &
  & 2*KF_LT, R_NDGNH, TDZAA, &
  & 1.0_JPRBT, &
  & ZINP, 2*KF_LT, TDZAA,&
  & ZAA, R_NDGNH, TDZAA, &
  & 0._JPRBT, &
  & ZOUTA, 2*KF_LT, R_NDGNH, &
  & D_NUMP)
#ifdef OMPGPU
!$OMP END TARGET DATA
#endif
#ifdef ACCGPU
!$ACC END HOST_DATA
#endif

! 2. +++++++++++++ symmetric
!IF KM=0 and NSMAX is 6:
!    IS=2
!    DO=1,4
!       PIA_2=2+1+(0..3)*2 ... 3+(0..3)*2 ... 3,5,7,9
!IF KM=0 and NSMAX is 7:
!    IS=1
!    DO=1,5
!       PIA_2=1+1+(1..5-1)*2 ...2+(0..4)*2 .... 2,4,6,8,10


#ifdef OMPGPU
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) PRIVATE(KM,ISKIP,ILS,IS) &
!$OMP&      SHARED(D_NUMP,R_NSMAX,KFC,KSTA,D_MYMS,IZBS,ITDZBS,ILDZBS,ZIA) DEFAULT(NONE)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM,IS,J) DEFAULT(NONE) &
!$ACC&    FIRSTPRIVATE(KF_LT,TDZAS)
#endif
DO KMLOC=1,D_NUMP
  DO JK=1,2*KF_LT
    KM =  D_MYMS(KMLOC)
    IS = 1+MOD(R_NSMAX-KM+1,2)
    IF(KM /= 0)THEN
      !$ACC LOOP SEQ
      DO J=1,(R_NSMAX-KM+3)/2
       ZINP(JK+(J-1+(KMLOC-1)*TDZAS)*2*KF_LT)=PIA(JK,IS+1+(J-1)*2,KMLOC)
      ENDDO
    ELSEIF (MOD((JK+1),2) .EQ. 0) THEN
      !$ACC LOOP SEQ
      DO J=1,(R_NSMAX+3)/2
        ZINP((JK-1)/2+1+(J-1+(KMLOC-1)*TDZAS)*2*KF_LT)=PIA(JK,IS+1+(J-1)*2,KMLOC)
      ENDDO
    ENDIF
  ENDDO
ENDDO

!C=A*B =>
! C^T=B^T*A^T

#ifdef OMPGPU
!$OMP TARGET DATA USE_DEVICE_PTR(ZAS,IZBS,IZCST)
#endif
#ifdef ACCGPU
!$ACC HOST_DATA USE_DEVICE(ZAS,ZINP,ZOUTS)
#endif
CALL HIP_GEMM_BATCHED( &
 & 'N', 'T', &
 & 2*KF_LT, R_NDGNH, TDZAS, &
 & 1.0_JPRBT, &
 & ZINP, 2*KF_LT, TDZAS, &
 & ZAS, R_NDGNH, TDZAS, &
 & 0._JPRBT, &
 & ZOUTS, 2*KF_LT, R_NDGNH, &
 & D_NUMP)
#ifdef OMPGPU
!$OMP END TARGET DATA
#endif
#ifdef ACCGPU
!$ACC END HOST_DATA
#endif


#ifdef OMPGPU
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) PRIVATE(KM,ISKIP,ISL) &
!$OMP&      SHARED(D_NUMP,R_NDGNH,KFC,D_MYMS,G_NDGLU,PSOA1,IZCST,ITDZCS,ILDZCS) DEFAULT(NONE)
#endif
#ifdef ACCGPU
!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM,ISL,IGLS,OFFSET1,OFFSET2,ISTAS,ZAOA,ZSOA) DEFAULT(NONE) &
!$ACC&    FIRSTPRIVATE(KF_LT) PRESENT(R_NDGL)
#endif
DO KMLOC=1,D_NUMP
  DO JK=1,2*KF_LT
    KM = D_MYMS(KMLOC)
    ISL = R_NDGNH-G_NDGLU(KM)+1
    !$ACC LOOP SEQ
    DO JGL=ISL,R_NDGNH
      IGLS = R_NDGL+1-JGL
      OFFSET1 = (D_NSTAGT1B(D_NPROCL(JGL)) + D_NPNTGTB1(KMLOC,JGL))*2*KF_LT
      OFFSET2 = (D_NSTAGT1B(D_NPROCL(IGLS)) + D_NPNTGTB1(KMLOC,IGLS))*2*KF_LT

      IF(KM /= 0) THEN
        ZSOA = ZOUTS(JK+(JGL-ISL+(KMLOC-1)*R_NDGNH)*2*KF_LT)
        ZAOA = ZOUTA(JK+(JGL-ISL+(KMLOC-1)*R_NDGNH)*2*KF_LT)
      ELSEIF (MOD((JK-1),2) .EQ. 0) THEN
        ZSOA = ZOUTS((JK-1)/2+1+(JGL-ISL+(KMLOC-1)*R_NDGNH)*2*KF_LT)
        ZAOA = ZOUTA((JK-1)/2+1+(JGL-ISL+(KMLOC-1)*R_NDGNH)*2*KF_LT)
       ENDIF

      FOUBUF_IN(OFFSET1+JK) = ZAOA+ZSOA
      FOUBUF_IN(OFFSET2+JK) = ZSOA-ZAOA
    ENDDO
  ENDDO
ENDDO
!$ACC END DATA

DEALLOCATE(ZINP)
DEALLOCATE(ZOUTS)
DEALLOCATE(ZOUTA)

IF (LSYNC_TRANS) THEN
  CALL MPL_BARRIER(CDSTRING='LEINV BARRIER')
ENDIF
CALL GSTATS(453,1)


IF (LHOOK) CALL DR_HOOK('LE_DGEMM',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE LEINV
END MODULE LEINV_MOD
