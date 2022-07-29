! (C) Copyright 2000- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TRMTOL_PACK_UNPACK
  USE ALLOCATOR_MOD
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: TRMTOL_PACK, TRMTOL_PACK_HANDLE, PREPARE_TRMTOL_PACK
  PUBLIC :: TRMTOL_UNPACK, TRMTOL_UNPACK_HANDLE, PREPARE_TRMTOL_UNPACK

  TYPE TRMTOL_PACK_HANDLE
    TYPE(ALLOCATION_RESERVATION_HANDLE) :: HFOUBUF_IN
  END TYPE
  TYPE TRMTOL_UNPACK_HANDLE
    TYPE(ALLOCATION_RESERVATION_HANDLE) :: HREEL
  END TYPE

CONTAINS
  FUNCTION PREPARE_TRMTOL_PACK(ALLOCATOR,KF_LEG) RESULT(HTRMTOL_PACK)
    USE PARKIND_ECTRANS, ONLY: JPIM, JPRBT, JPRD
    USE TPM_DISTR, ONLY: D
    USE TPM_DIM, ONLY: R
    USE ISO_C_BINDING
    USE LEINV_MOD
    USE ALLOCATOR_MOD

    IMPLICIT NONE

    TYPE(BUFFERED_ALLOCATOR), INTENT(INOUT) :: ALLOCATOR
    INTEGER(KIND=JPIM), INTENT(IN) :: KF_LEG

    TYPE(TRMTOL_PACK_HANDLE) :: HTRMTOL_PACK

    INTEGER(KIND=C_SIZE_T) :: IALLOC_SZ

    REAL(KIND=JPRBT) :: ZPRBT_DUMMY

    IALLOC_SZ = D%NLENGT1B*2*KF_LEG*SIZEOF(ZPRBT_DUMMY)
    HTRMTOL_PACK%HFOUBUF_IN = RESERVE(ALLOCATOR, IALLOC_SZ)
  END FUNCTION
  SUBROUTINE TRMTOL_PACK(ALLOCATOR,HTRMTOL_PACK,ZOUTS,ZOUTA,ZOUTS0,ZOUTA0,FOUBUF_IN,KF_LEG)

    !**** *TRMTOL_PACK* - Packing buffer for TRMTOL

    !     Purpose.
    !     --------
    !        Packs data from LTINV outputs into FOUBUF for conversion to fourier space

    !**   Interface.
    !     ----------
    !        CALL TRMTOL_PACK(...)

    !        Explicit arguments :  ZOUTS - symmetric data
    !        --------------------  ZOUTA - asymmetric data
    !                              ZOUTS0 - symmetric data for KMLOC0
    !                              ZOUTA0 - asymmetric data for KMLOC0
    !                              FOUBUF_IN - output towards TRMTOL
    !                              KF_LEG - number of fields (we have 2XKF_LEG because complex)

    !        Implicit arguments :  None.
    !        --------------------

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

    USE PARKIND_ECTRANS  ,ONLY : JPIM,JPRB,JPRBT,JPRD
    USE YOMHOOK, ONLY : LHOOK,DR_HOOK, JPHOOK
    USE TPM_DIM, ONLY : R, R_NDGNH,R_NDGL
    USE TPM_GEOMETRY,ONLY : G,G_NDGLU
    USE TPM_DISTR, ONLY : D,D_NUMP,D_MYMS,D_NPNTGTB1
    USE LEINV_MOD, ONLY: LEINV_STRIDES
    USE TPM_TRANS, ONLY: LEINV_CONF

    IMPLICIT NONE


    !     DUMMY ARGUMENTS
    TYPE(BUFFERED_ALLOCATOR), INTENT(IN) :: ALLOCATOR
    TYPE(TRMTOL_PACK_HANDLE), INTENT(IN) :: HTRMTOL_PACK
    REAL(KIND=JPRB), INTENT(OUT), POINTER :: FOUBUF_IN(:)
    REAL(KIND=JPRBT), INTENT(IN) :: ZOUTS(:), ZOUTA(:)
    REAL(KIND=JPRD), INTENT(IN) :: ZOUTS0(:), ZOUTA0(:)
    INTEGER(KIND=JPIM), INTENT(IN)  :: KF_LEG

    !     LOCAL
    REAL(KIND=JPRBT) :: ZAOA, ZSOA

    INTEGER(KIND=JPIM) :: KMLOC, KM, ISL, JGL, JK, IGLS, OFFSET1, OFFSET2
    INTEGER(KIND=JPIM)  :: IOUT_STRIDES0, IOUT_SIZE
    INTEGER(KIND=JPIM)  :: IOUT0_STRIDES0, IOUT0_SIZE

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('TRMTOL_PACK',0,ZHOOK_HANDLE)

    CALL ASSIGN_PTR(FOUBUF_IN, GET_ALLOCATION(ALLOCATOR, HTRMTOL_PACK%HFOUBUF_IN),&
      & 1_C_SIZE_T, D%NLENGT1B*2*KF_LEG*SIZEOF(FOUBUF_IN(1)))

    CALL LEINV_STRIDES(KF_LEG,IOUT_STRIDES0=IOUT_STRIDES0,IOUT_SIZE=IOUT_SIZE,&
                       IOUT0_STRIDES0=IOUT0_STRIDES0,IOUT0_SIZE=IOUT0_SIZE)

    !$ACC DATA PRESENT(D,D_MYMS,D_NPNTGTB1,G,G_NDGLU) &
    !$ACC&     PRESENT(ZOUTS,ZOUTA,ZOUTS0,ZOUTA0,FOUBUF_IN,LEINV_CONF)

    !$ACC PARALLEL LOOP COLLAPSE(3) DEFAULT(NONE) PRIVATE(KM,ISL,IGLS,OFFSET1,OFFSET2,ZAOA,ZSOA) ASYNC(1)
    DO KMLOC=1,D_NUMP
      DO JGL=1,R_NDGNH
        DO JK=1,2*KF_LEG
          KM = D_MYMS(KMLOC)
          ISL = R_NDGNH-G_NDGLU(KM)+1
          IF (JGL >= ISL) THEN
            !(DO JGL=ISL,R_NDGNH)
            IGLS = R_NDGL+1-JGL
            OFFSET1 = D_NPNTGTB1(KMLOC,JGL )*2*KF_LEG
            OFFSET2 = D_NPNTGTB1(KMLOC,IGLS)*2*KF_LEG

            IF(KM /= 0) THEN
              ZSOA = ZOUTS(JK+(JGL-ISL)*IOUT_STRIDES0+LEINV_CONF%OFFSETS_N(KMLOC)*IOUT_STRIDES0)
              ZAOA = ZOUTA(JK+(JGL-ISL)*IOUT_STRIDES0+LEINV_CONF%OFFSETS_N(KMLOC)*IOUT_STRIDES0)
            ELSEIF (MOD((JK-1),2) .EQ. 0) THEN
              ZSOA = ZOUTS0((JK-1)/2+1+(JGL-1)*IOUT0_STRIDES0)
              ZAOA = ZOUTA0((JK-1)/2+1+(JGL-1)*IOUT0_STRIDES0)
            ELSE
              ! Imaginary values of KM=0 is zero, though I don't think we care
              ZSOA = 0_JPRBT
              ZAOA = 0_JPRBT
            ENDIF

            FOUBUF_IN(OFFSET1+JK) = ZAOA+ZSOA
            FOUBUF_IN(OFFSET2+JK) = ZSOA-ZAOA
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    !$ACC WAIT(1)

    !$ACC END DATA

    IF (LHOOK) CALL DR_HOOK('TRMTOL_PACK',1,ZHOOK_HANDLE)

  END SUBROUTINE TRMTOL_PACK

  FUNCTION PREPARE_TRMTOL_UNPACK(ALLOCATOR,KF_FS) RESULT(HTRMTOL_UNPACK)
    USE PARKIND_ECTRANS, ONLY: JPIM, JPRBT
    USE TPM_DISTR, ONLY: D

    IMPLICIT NONE

    TYPE(BUFFERED_ALLOCATOR), INTENT(INOUT) :: ALLOCATOR
    INTEGER(KIND=JPIM) :: KF_FS

    TYPE(TRMTOL_UNPACK_HANDLE) :: HTRMTOL_UNPACK

    REAL(KIND=JPRBT) :: DUMMY

    HTRMTOL_UNPACK%HREEL = RESERVE(ALLOCATOR, D%NLENGTF*KF_FS*SIZEOF(DUMMY))

  END FUNCTION
SUBROUTINE TRMTOL_UNPACK(ALLOCATOR,HTRMTOL_UNPACK,FOUBUF,PREEL_COMPLEX,KF_CURRENT,KF_TOTAL)

!**** *TRMTOL_UNPACK* - Copy fourier data from buffer to local array

!     Purpose.
!     --------
!        Routine for copying fourier data from buffer to local array

!**   Interface.
!     ----------
!     CALL TRMTOL_UNPACK(...)

!     Explicit arguments :  PREEL_COMPLEX - local fourier/GP array
!     --------------------  KF_CURRENT - number of fields that are read (from Legendre space)
!                           KF_TOTAL - total fields in PREEL ("stride")
!
!     Externals.  None.
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 2000-04-01

!     ------------------------------------------------------------------

USE PARKIND_ECTRANS  ,ONLY : JPIM,JPRBT
USE TPM_DISTR       ,ONLY : D,MYSETW,MYPROC, NPROC, D_NSTAGTF, D_NPNTGTB0,D_NPTRLS
USE TPM_GEOMETRY    ,ONLY : G,G_NMEN,G_NLOEN,G_NLOEN_MAX
!

IMPLICIT NONE

REAL(KIND=JPRBT), INTENT(IN) :: FOUBUF(:)
REAL(KIND=JPRBT), INTENT(OUT), POINTER :: PREEL_COMPLEX(:)
INTEGER(KIND=JPIM),INTENT(IN) :: KF_CURRENT, KF_TOTAL
TYPE(BUFFERED_ALLOCATOR), INTENT(IN) :: ALLOCATOR
TYPE(TRMTOL_UNPACK_HANDLE), INTENT(IN) :: HTRMTOL_UNPACK

INTEGER(KIND=JPIM) :: JM,JF,IGLG,ISTA,OFFSET_VAR,IOFF_LAT,KGL
REAL(KIND=JPRBT) :: RET_REAL, RET_COMPLEX

CALL ASSIGN_PTR(PREEL_COMPLEX, GET_ALLOCATION(ALLOCATOR, HTRMTOL_UNPACK%HREEL),&
    & 1_C_SIZE_T, KF_TOTAL*D%NLENGTF*SIZEOF(PREEL_COMPLEX(1)))

!$ACC DATA PRESENT(D,G_NLOEN,G_NMEN,D_NPNTGTB0,FOUBUF,PREEL_COMPLEX,D_NSTAGTF) ASYNC(1)

OFFSET_VAR=D_NPTRLS(MYSETW)
!$ACC PARALLEL LOOP PRIVATE(IGLG,IOFF_LAT,ISTA,RET_REAL,RET_COMPLEX) DEFAULT(NONE) &
!$ACC& ASYNC(1) TILE(32,16,1)
DO KGL=1,D%NDGL_FS
  DO JF=1,KF_CURRENT
    DO JM=0,G_NLOEN_MAX/2
      IGLG = OFFSET_VAR+KGL-1

      ! FFT transforms NLON real values to floor(NLON/2)+1 complex numbers. Hence we have
      ! to fill those floor(NLON/2)+1 values.
      ! Truncation happens starting at G_NMEN+1. Hence, we zero-fill those values.
      IF (JM <= G_NLOEN(IGLG)/2) THEN
        RET_REAL = 0.0_JPRBT
        RET_COMPLEX = 0.0_JPRBT
        IF (JM <= G_NMEN(IGLG)) THEN
          ISTA  = D_NPNTGTB0(JM,KGL)*KF_CURRENT*2

          RET_REAL    = FOUBUF(ISTA+2*JF-1)
          RET_COMPLEX = FOUBUF(ISTA+2*JF  )
        ENDIF
        IOFF_LAT = KF_TOTAL*D_NSTAGTF(KGL)+(JF-1)*(D_NSTAGTF(KGL+1)-D_NSTAGTF(KGL))
        PREEL_COMPLEX(IOFF_LAT+2*JM+1) = RET_REAL
        PREEL_COMPLEX(IOFF_LAT+2*JM+2) = RET_COMPLEX
      ENDIF
    ENDDO
  ENDDO
ENDDO
!$ACC END DATA

!$ACC WAIT(1)

END SUBROUTINE TRMTOL_UNPACK
END MODULE TRMTOL_PACK_UNPACK

