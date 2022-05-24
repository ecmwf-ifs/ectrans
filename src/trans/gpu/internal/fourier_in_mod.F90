! (C) Copyright 2000- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FOURIER_IN_MOD
  USE ALLOCATOR_MOD
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: FOURIER_IN, FOURIER_IN_HANDLE, PREPARE_FOURIER_IN

  TYPE FOURIER_IN_HANDLE
    TYPE(ALLOCATION_RESERVATION_HANDLE) :: HREEL
  END TYPE

CONTAINS
  FUNCTION PREPARE_FOURIER_IN(ALLOCATOR,KF_FS) RESULT(HFOURIER_IN)
    USE PARKIND_ECTRANS, ONLY: JPIM, JPRBT
    USE TPM_DISTR, ONLY: D

    IMPLICIT NONE

    TYPE(BUFFERED_ALLOCATOR), INTENT(INOUT) :: ALLOCATOR
    INTEGER(KIND=JPIM) :: KF_FS

    TYPE(FOURIER_IN_HANDLE) :: HFOURIER_IN

    REAL(KIND=JPRBT) :: DUMMY

    HFOURIER_IN%HREEL = RESERVE(ALLOCATOR, D%NLENGTF*KF_FS*SIZEOF(DUMMY))

  END FUNCTION
SUBROUTINE FOURIER_IN(ALLOCATOR,HFOURIER_IN,FOUBUF,PREEL_COMPLEX,KF_CURRENT,KF_TOTAL)

!**** *FOURIER_IN* - Copy fourier data from buffer to local array

!     Purpose.
!     --------
!        Routine for copying fourier data from buffer to local array

!**   Interface.
!     ----------
!     CALL FOURIER_IN(...)

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
TYPE(FOURIER_IN_HANDLE), INTENT(IN) :: HFOURIER_IN

INTEGER(KIND=JPIM) :: JM,JF,IGLG,ISTA,OFFSET_VAR,IOFF_LAT,KGL
REAL(KIND=JPRBT) :: RET_REAL, RET_COMPLEX

CALL ASSIGN_PTR(PREEL_COMPLEX, GET_ALLOCATION(ALLOCATOR, HFOURIER_IN%HREEL),&
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

END SUBROUTINE FOURIER_IN
END MODULE FOURIER_IN_MOD

