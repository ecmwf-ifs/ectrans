! (C) Copyright 2022- NVIDIA.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
#define ALIGN(I, A) (((I)+(A)-1)/(A)*(A))
MODULE BUFFERED_ALLOCATOR_MOD

  USE EC_PARKIND,            ONLY: JPIM
  USE ABORT_TRANS_MOD,       ONLY: ABORT_TRANS
  USE ISO_C_BINDING,         ONLY: C_INT8_T, C_SIZE_T, C_LOC, C_F_POINTER
  USE GROWING_ALLOCATOR_MOD, ONLY: GROWING_ALLOCATION_TYPE
  USE OPENACC,               ONLY: ACC_ASYNC_SYNC

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: BUFFERED_ALLOCATOR, ALLOCATION_RESERVATION_HANDLE, RESERVE, ASSIGN_PTR, GET_ALLOCATION
  PUBLIC :: MAKE_BUFFERED_ALLOCATOR, INSTANTIATE_ALLOCATOR

  ! The buffered allocator uses double buffering. The idea is that the allocator
  ! iterates through its two buffers, and each allocate returns one or the other
  ! buffer. It is a two-step allocator - it expects you to create reservation
  ! handles first for all allocations. Then the allocator is instantiated (i.e.
  ! the buffers are actually allocated). Instantiation will do an allocation
  ! that is large enough two hold all consecutive allocations. Other allocations
  ! might be overwritten (like you can't access the allocation done two steps
  ! before).
  ! After instantiation, you can retrieve your buffers by passing the allocator
  ! and the handles to GET_ALLOCATION. Also, we provide helper function
  ! ASSIGN_PTR, because an allocation is often split among several "sub-buffers",
  ! so you can for example assign the first half of an allocation to one
  ! buffer, while the second half to another buffer.
  ! If you see "Logical errors" that usually means you try to retrieve a buffer
  ! that is not within the reserved allocation size. This might be a valid
  ! region in the sense that it is physically allocated, but it might be part of
  ! the double buffer.


  INTEGER(KIND=JPIM), PARAMETER :: NBUF = 2
  TYPE BUFFERED_ALLOCATOR
    INTEGER(KIND=C_SIZE_T) :: BUFR_SZ(0:NBUF-1)
    INTEGER(KIND=JPIM) :: NEXT_BUF
    TYPE(GROWING_ALLOCATION_TYPE), POINTER :: PTR
  END TYPE
  TYPE ALLOCATION_RESERVATION_HANDLE
    INTEGER(KIND=C_SIZE_T) :: SZ
    INTEGER(KIND=JPIM) :: BUF
  END TYPE

  INTERFACE ASSIGN_PTR
    MODULE PROCEDURE ASSIGN_PTR_FLOAT, ASSIGN_PTR_DOUBLE
  END INTERFACE

CONTAINS

  ! TODO This is not perfect yet. We will over-allocate up to 2X in theory.
  ! It would be better to always keep the previous allocation size and then
  ! have one allocation sitting at the the top, and the double-buffer at
  ! the bottom of the allocation.

  FUNCTION MAKE_BUFFERED_ALLOCATOR()
    IMPLICIT NONE
    TYPE(BUFFERED_ALLOCATOR) :: MAKE_BUFFERED_ALLOCATOR

    MAKE_BUFFERED_ALLOCATOR%BUFR_SZ(:) = 0
    MAKE_BUFFERED_ALLOCATOR%NEXT_BUF = 0
  END FUNCTION MAKE_BUFFERED_ALLOCATOR

  FUNCTION RESERVE(ALLOCATOR, SZ)
    IMPLICIT NONE
    TYPE(BUFFERED_ALLOCATOR), INTENT(INOUT) :: ALLOCATOR
    INTEGER(KIND=C_SIZE_T), INTENT(IN) :: SZ

    TYPE(ALLOCATION_RESERVATION_HANDLE) :: RESERVE

    ALLOCATOR%BUFR_SZ(ALLOCATOR%NEXT_BUF) = MAX(ALLOCATOR%BUFR_SZ(ALLOCATOR%NEXT_BUF),SZ)
    RESERVE%BUF = ALLOCATOR%NEXT_BUF
    RESERVE%SZ = SZ

    ALLOCATOR%NEXT_BUF = MOD(ALLOCATOR%NEXT_BUF+1,NBUF)
  END FUNCTION RESERVE

  SUBROUTINE INSTANTIATE_ALLOCATOR(ALLOCATOR, GROWING_ALLOCATION)
    USE GROWING_ALLOCATOR_MOD, ONLY: REALLOCATE_GROWING_ALLOCATION
    IMPLICIT NONE
    TYPE(BUFFERED_ALLOCATOR), INTENT(INOUT) :: ALLOCATOR
    !!TYPE(GROWING_ALLOCATION_TYPE), INTENT(IN), POINTER :: GROWING_ALLOCATION
    TYPE(GROWING_ALLOCATION_TYPE), target, INTENT(INout) :: GROWING_ALLOCATION
    INTEGER :: I

    DO I = 0, NBUF-1
      ALLOCATOR%BUFR_SZ(I) = ALIGN(ALLOCATOR%BUFR_SZ(I),128)
    ENDDO
    ALLOCATOR%PTR => GROWING_ALLOCATION

    CALL REALLOCATE_GROWING_ALLOCATION(GROWING_ALLOCATION, SUM(ALLOCATOR%BUFR_SZ))
  END SUBROUTINE

  FUNCTION GET_ALLOCATION(ALLOCATOR, RESERVATION)
    IMPLICIT NONE
    TYPE(BUFFERED_ALLOCATOR), INTENT(IN) :: ALLOCATOR
    TYPE(ALLOCATION_RESERVATION_HANDLE), INTENT(IN) :: RESERVATION

    INTEGER(KIND=C_INT8_T), POINTER :: GET_ALLOCATION(:)

    IF (RESERVATION%SZ > ALLOCATOR%BUFR_SZ(RESERVATION%BUF)) THEN
      CALL ABORT_TRANS( "Logical Error in GET_ALLOCATION")
    ENDIF
    IF (RESERVATION%BUF == 0) THEN
      GET_ALLOCATION(1:) => ALLOCATOR%PTR%PTR(1:RESERVATION%SZ)
    ELSE
      GET_ALLOCATION(1:) => ALLOCATOR%PTR%PTR(SUM(ALLOCATOR%BUFR_SZ(0:RESERVATION%BUF-1))+1: &
                                              SUM(ALLOCATOR%BUFR_SZ(0:RESERVATION%BUF-1))+RESERVATION%SZ)
    ENDIF
  END FUNCTION GET_ALLOCATION

  SUBROUTINE ASSIGN_PTR_FLOAT(DST, SRC, START_IN_BYTES, LENGTH_IN_BYTES, SET_VALUE, SET_STREAM)
    USE ISO_C_BINDING, ONLY: C_FLOAT
    IMPLICIT NONE
    INTEGER(KIND=C_INT8_T), POINTER, INTENT(INOUT) :: SRC(:)
    REAL(KIND=C_FLOAT), POINTER, INTENT(OUT) :: DST(:)
    LOGICAL, INTENT(IN), OPTIONAL :: SET_VALUE
    INTEGER(KIND=4), INTENT(IN), OPTIONAL :: SET_STREAM
    LOGICAL :: SET_VALUE_EFF
    INTEGER(KIND=4) :: SET_STREAM_EFF
    INTEGER(KIND=C_SIZE_T) :: START_IN_BYTES, LENGTH_IN_BYTES
    IF (START_IN_BYTES + LENGTH_IN_BYTES - 1 > SIZE(SRC, KIND=C_SIZE_T)) THEN
      CALL ABORT_TRANS("Logical Error in ASSIGN_PTR - OOB assignment")
    ENDIF
    IF (PRESENT(SET_VALUE)) THEN
      SET_VALUE_EFF = SET_VALUE
    ELSE
      SET_VALUE_EFF = .FALSE.
    ENDIF
    IF (PRESENT(SET_STREAM)) THEN
        SET_STREAM_EFF = SET_STREAM
    ELSE
        SET_STREAM_EFF = ACC_ASYNC_SYNC
    ENDIF
    IF (SET_VALUE_EFF .AND. LENGTH_IN_BYTES > 0) THEN
      ! This option is turned off by default, but for experimentation we can turn it on. This is
      ! setting all bits to 1 (meaning NaN in floating point)
      !$ACC KERNELS PRESENT(SRC) ASYNC(SET_STREAM_EFF)
      SRC(START_IN_BYTES:START_IN_BYTES+LENGTH_IN_BYTES-1) = -1
      !$ACC END KERNELS!! LOOP
    ENDIF
    CALL C_F_POINTER(C_LOC(SRC(START_IN_BYTES:START_IN_BYTES+LENGTH_IN_BYTES-1)), DST, &
        & [STORAGE_SIZE(SRC(START_IN_BYTES:START_IN_BYTES+LENGTH_IN_BYTES-1))/STORAGE_SIZE(DST(0))])
  END SUBROUTINE ASSIGN_PTR_FLOAT
  SUBROUTINE ASSIGN_PTR_DOUBLE(DST, SRC, START_IN_BYTES, LENGTH_IN_BYTES, SET_VALUE, SET_STREAM)
    USE ISO_C_BINDING, ONLY: C_DOUBLE
    IMPLICIT NONE
    INTEGER(KIND=C_INT8_T), POINTER, INTENT(INOUT) :: SRC(:)
    REAL(KIND=C_DOUBLE), POINTER, INTENT(OUT) :: DST(:)
    LOGICAL, INTENT(IN), OPTIONAL :: SET_VALUE
    INTEGER(KIND=4), INTENT(IN), OPTIONAL :: SET_STREAM
    LOGICAL :: SET_VALUE_EFF
    INTEGER(KIND=4) :: SET_STREAM_EFF
    INTEGER(KIND=C_SIZE_T) :: START_IN_BYTES, LENGTH_IN_BYTES
    IF (START_IN_BYTES + LENGTH_IN_BYTES - 1 > SIZE(SRC, KIND=C_SIZE_T)) THEN
      CALL ABORT_TRANS("Logical Error in ASSIGN_PTR - OOB assignment")
    ENDIF
    IF (PRESENT(SET_VALUE)) THEN
      SET_VALUE_EFF = SET_VALUE
    ELSE
      SET_VALUE_EFF = .FALSE.
    ENDIF
    IF (PRESENT(SET_STREAM)) THEN
        SET_STREAM_EFF = SET_STREAM
    ELSE
        SET_STREAM_EFF = ACC_ASYNC_SYNC
    ENDIF
    IF (SET_VALUE_EFF .AND. LENGTH_IN_BYTES > 0) THEN
      ! This option is turned off by default, but for experimentation we can turn it on. This is
      ! setting all bits to 1 (meaning NaN in floating point)
      !$ACC KERNELS PRESENT(SRC) ASYNC(SET_STREAM_EFF)
      SRC(START_IN_BYTES:START_IN_BYTES+LENGTH_IN_BYTES-1) = -1
      !$ACC END KERNELS!! LOOP
    ENDIF
    CALL C_F_POINTER(C_LOC(SRC(START_IN_BYTES:START_IN_BYTES+LENGTH_IN_BYTES-1)), DST, &
        & [STORAGE_SIZE(SRC(START_IN_BYTES:START_IN_BYTES+LENGTH_IN_BYTES-1))/STORAGE_SIZE(DST(0))])
  END SUBROUTINE ASSIGN_PTR_DOUBLE
END MODULE
