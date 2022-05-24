! (C) Copyright 2000- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FOURIER_OUT_MOD
  USE ALLOCATOR_MOD
  IMPLICIT NONE

  TYPE FOURIER_OUT_HANDLE
    TYPE(ALLOCATION_RESERVATION_HANDLE) :: HFOUBUF_IN
  END TYPE
CONTAINS
  FUNCTION PREPARE_FOURIER_OUT(ALLOCATOR, KF_FS) RESULT(HFOURIER_OUT)
    USE PARKIND_ECTRANS, ONLY: JPIM, JPRBT
    USE TPM_DISTR, ONLY: D

    IMPLICIT NONE

    TYPE(BUFFERED_ALLOCATOR), INTENT(INOUT) :: ALLOCATOR
    INTEGER(KIND=JPIM), INTENT(IN) :: KF_FS
    TYPE(FOURIER_OUT_HANDLE) :: HFOURIER_OUT

    REAL(KIND=JPRBT) :: DUMMY

    HFOURIER_OUT%HFOUBUF_IN = RESERVE(ALLOCATOR, D%NLENGT0B*KF_FS*2*SIZEOF(DUMMY))
  END FUNCTION
SUBROUTINE FOURIER_OUT(ALLOCATOR,HFOURIER_OUT,PREEL_COMPLEX,FOUBUF_IN,KF_FS)

!**** *FOURIER_OUT* - Copy fourier data from local array to buffer

!     Purpose.
!     --------
!        Routine for copying fourier data from local array to buffer

!**   Interface.
!     ----------
!     CALL FOURIER_OUT(...)

!     Explicit arguments :  PREEL - local fourier/GP array
!     --------------------  KF_FS - number of fields
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

USE ALLOCATOR_MOD
USE PARKIND_ECTRANS  ,ONLY : JPIM,JPRBT
USE TPM_DISTR       ,ONLY : D,MYSETW,MYPROC, NPROC, D_NSTAGTF, D_NPNTGTB0,D_NPTRLS
USE TPM_GEOMETRY    ,ONLY : G,G_NMEN,G_NLOEN
USE TPM_DIM, ONLY: R_NSMAX
USE ISO_C_BINDING
!

IMPLICIT NONE

REAL(KIND=JPRBT), INTENT(IN) :: PREEL_COMPLEX(:)
REAL(KIND=JPRBT), POINTER, INTENT(OUT) :: FOUBUF_IN(:)
INTEGER(KIND=JPIM),INTENT(IN) :: KF_FS
TYPE(BUFFERED_ALLOCATOR), INTENT(IN) :: ALLOCATOR
TYPE(FOURIER_OUT_HANDLE), INTENT(IN) :: HFOURIER_OUT

INTEGER(KIND=JPIM) :: JM,JF,IGLG,ISTA,OFFSET_VAR,IOFF_LAT,KGL

REAL(KIND=JPRBT)    :: SCAL

CALL ASSIGN_PTR(FOUBUF_IN, GET_ALLOCATION(ALLOCATOR, HFOURIER_OUT%HFOUBUF_IN),&
    & 1_C_SIZE_T, D%NLENGT0B*KF_FS*2*SIZEOF(FOUBUF_IN(1)))

!$ACC DATA PRESENT(D,G_NMEN,D_NPNTGTB0,FOUBUF_IN,PREEL_COMPLEX,D_NSTAGTF,G_NLOEN) ASYNC(1)

! scale results and move into next transformation buffer

OFFSET_VAR=D_NPTRLS(MYSETW)

!$ACC PARALLEL LOOP PRIVATE(IGLG,IOFF_LAT,ISTA,SCAL) DEFAULT(NONE) &
!$ACC& ASYNC(1) TILE(32,16,1)
DO KGL=1,D%NDGL_FS
  DO JM=0,R_NSMAX !(note that R_NSMAX <= G_NMEN(IGLG) for all IGLG)
    DO JF=1,KF_FS
      IGLG = OFFSET_VAR+KGL-1
      IF (JM <= G_NMEN(IGLG)) THEN
        IOFF_LAT = KF_FS*D_NSTAGTF(KGL)+(JF-1)*(D_NSTAGTF(KGL+1)-D_NSTAGTF(KGL))

        SCAL = 1._JPRBT/REAL(G_NLOEN(IGLG),JPRBT)
        ISTA  = D_NPNTGTB0(JM,KGL)*KF_FS*2

        FOUBUF_IN(ISTA+2*JF-1) = SCAL * PREEL_COMPLEX(IOFF_LAT+2*JM+1)
        FOUBUF_IN(ISTA+2*JF  ) = SCAL * PREEL_COMPLEX(IOFF_LAT+2*JM+2)
      ENDIF
    ENDDO
  ENDDO
ENDDO
!$ACC END DATA

!$ACC WAIT(1)

END SUBROUTINE FOURIER_OUT
END MODULE FOURIER_OUT_MOD

