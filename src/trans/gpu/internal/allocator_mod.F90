#define ALIGN(I, A) (((I)+(A)-1)/(A)*(A))
MODULE ALLOCATOR_MOD

  USE PARKIND_ECTRANS ,ONLY : JPIM
  USE ISO_C_BINDING, ONLY: C_INT8_T, C_SIZE_T

  IMPLICIT NONE

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


  TYPE BUFFERED_ALLOCATOR
    INTEGER(KIND=C_SIZE_T) :: BUFR_SZ(0:1)
    INTEGER(KIND=JPIM) :: NEXT_BUF
    INTEGER(C_INT8_T), POINTER :: PTR(:)
  END TYPE
  TYPE ALLOCATION_RESERVATION_HANDLE
    INTEGER(KIND=C_SIZE_T) :: SZ
    INTEGER(KIND=JPIM) :: BUF
  END TYPE

  INTERFACE ASSIGN_PTR
    SUBROUTINE ASSIGN_PTR_FLOAT(DST, SRC, START_IN_BYTES, LENGTH_IN_BYTES)
      USE ISO_C_BINDING
      INTEGER(KIND=C_INT8_T), POINTER, INTENT(IN) :: SRC(:)
      REAL(KIND=C_FLOAT), POINTER, INTENT(OUT) :: DST(:)
      INTEGER(KIND=C_SIZE_T) :: START_IN_BYTES, LENGTH_IN_BYTES
    END SUBROUTINE
    SUBROUTINE ASSIGN_PTR_DOUBLE(DST, SRC, START_IN_BYTES, LENGTH_IN_BYTES)
      USE ISO_C_BINDING
      INTEGER(KIND=C_INT8_T), POINTER, INTENT(IN) :: SRC(:)
      REAL(KIND=C_DOUBLE), POINTER, INTENT(OUT) :: DST(:)
      INTEGER(KIND=C_SIZE_T) :: START_IN_BYTES, LENGTH_IN_BYTES
    END SUBROUTINE
  END INTERFACE

CONTAINS

  ! TODO This is not perfect yet. We will over-allocate up to 2X in theory.
  ! It would be better to always keep the previous allocation size and then
  ! have one allocation sitting at the the top, and the double-buffer at
  ! the bottom of the allocation.

  FUNCTION MAKE_BUFFERED_ALLOCATOR()
    TYPE(BUFFERED_ALLOCATOR) :: MAKE_BUFFERED_ALLOCATOR

    MAKE_BUFFERED_ALLOCATOR%BUFR_SZ(:) = 0
    MAKE_BUFFERED_ALLOCATOR%NEXT_BUF = 0
  END FUNCTION

  FUNCTION RESERVE(ALLOCATOR, SZ)
    TYPE(BUFFERED_ALLOCATOR), INTENT(INOUT) :: ALLOCATOR
    INTEGER(KIND=C_SIZE_T), INTENT(IN) :: SZ

    TYPE(ALLOCATION_RESERVATION_HANDLE) :: RESERVE

    ALLOCATOR%BUFR_SZ(ALLOCATOR%NEXT_BUF) = MAX(ALLOCATOR%BUFR_SZ(ALLOCATOR%NEXT_BUF),SZ)
    RESERVE%BUF = ALLOCATOR%NEXT_BUF
    RESERVE%SZ = SZ

    ALLOCATOR%NEXT_BUF = 1-ALLOCATOR%NEXT_BUF
  END FUNCTION

  SUBROUTINE INSTANTIATE_ALLOCATOR(ALLOCATOR, OLD_PTR)
    TYPE(BUFFERED_ALLOCATOR), INTENT(INOUT) :: ALLOCATOR
    INTEGER(C_INT8_T), OPTIONAL, INTENT(INOUT), POINTER :: OLD_PTR(:)

    ALLOCATOR%BUFR_SZ(1) = ALIGN(ALLOCATOR%BUFR_SZ(1),128)
    ALLOCATOR%BUFR_SZ(2) = ALIGN(ALLOCATOR%BUFR_SZ(2),128)

    IF (ASSOCIATED(OLD_PTR)) THEN
      IF (SIZEOF(OLD_PTR) < SUM(ALLOCATOR%BUFR_SZ) ) THEN
        !$ACC EXIT DATA DELETE(OLD_PTR) IF(PRESENT(OLD_PTR))
        DEALLOCATE(OLD_PTR)
        NULLIFY(OLD_PTR)
        ALLOCATE(ALLOCATOR%PTR(1:SUM(ALLOCATOR%BUFR_SZ)))
        !$ACC ENTER DATA CREATE(ALLOCATOR%PTR)

        OLD_PTR => ALLOCATOR%PTR
      ELSE
        ALLOCATOR%PTR(1:) => OLD_PTR(1:)
      ENDIF
    ELSE
      ALLOCATE(ALLOCATOR%PTR(1:SUM(ALLOCATOR%BUFR_SZ)))
      !$ACC ENTER DATA CREATE(ALLOCATOR%PTR)
      OLD_PTR => ALLOCATOR%PTR
    ENDIF
  END SUBROUTINE

  FUNCTION GET_ALLOCATION(ALLOCATOR, RESERVATION)
    TYPE(BUFFERED_ALLOCATOR), INTENT(IN) :: ALLOCATOR
    TYPE(ALLOCATION_RESERVATION_HANDLE), INTENT(IN) :: RESERVATION

    INTEGER(KIND=C_INT8_T), POINTER :: GET_ALLOCATION(:)

    IF (RESERVATION%SZ > ALLOCATOR%BUFR_SZ(RESERVATION%BUF)) THEN
      PRINT *, "Logical Error in GET_ALLOCATOIN"
      STOP 4
    ENDIF
    IF (RESERVATION%BUF == 0) THEN
      GET_ALLOCATION(1:) => ALLOCATOR%PTR(1:RESERVATION%SZ)
    ELSE
      GET_ALLOCATION(1:) => ALLOCATOR%PTR(ALLOCATOR%BUFR_SZ(0)+1: &
                                          ALLOCATOR%BUFR_SZ(0)+RESERVATION%SZ)
    ENDIF
  END FUNCTION

  SUBROUTINE ASSIGN_PTR_FLOAT(DST, SRC, START_IN_BYTES, LENGTH_IN_BYTES)
    USE ISO_C_BINDING
    INTEGER(KIND=C_INT8_T), POINTER, INTENT(IN) :: SRC(:)
    REAL(KIND=C_FLOAT), POINTER, INTENT(OUT) :: DST(:)
    INTEGER(KIND=C_SIZE_T) :: START_IN_BYTES, LENGTH_IN_BYTES
    IF (START_IN_BYTES + LENGTH_IN_BYTES - 1 > SIZE(SRC, KIND=C_SIZE_T)) THEN
      PRINT *, "Logical Error in ASSIGN_PTR - OOB assignment"
      STOP 4
    ENDIF
    CALL C_F_POINTER(C_LOC(SRC(START_IN_BYTES:START_IN_BYTES+LENGTH_IN_BYTES)), DST, &
        & [SIZEOF(SRC(START_IN_BYTES:START_IN_BYTES+LENGTH_IN_BYTES))/SIZEOF(DST(0))])
  END SUBROUTINE
  SUBROUTINE ASSIGN_PTR_DOUBLE(DST, SRC, START_IN_BYTES, LENGTH_IN_BYTES)
    USE ISO_C_BINDING
    INTEGER(KIND=C_INT8_T), POINTER, INTENT(IN) :: SRC(:)
    REAL(KIND=C_DOUBLE), POINTER, INTENT(OUT) :: DST(:)
    INTEGER(KIND=C_SIZE_T) :: START_IN_BYTES, LENGTH_IN_BYTES
    IF (START_IN_BYTES + LENGTH_IN_BYTES - 1 > SIZE(SRC, KIND=C_SIZE_T)) THEN
      PRINT *, "Logical Error in ASSIGN_PTR - OOB assignment"
      STOP 4
    ENDIF
    CALL C_F_POINTER(C_LOC(SRC(START_IN_BYTES:START_IN_BYTES+LENGTH_IN_BYTES)), DST, &
        & [SIZEOF(SRC(START_IN_BYTES:START_IN_BYTES+LENGTH_IN_BYTES))/SIZEOF(DST(0))])
  END SUBROUTINE
END MODULE
