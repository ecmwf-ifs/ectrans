!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                               !!
  !!                    FortsRaw                   !!
  !!                                               !!
  !!      Copyright (c) 2019, Thomas Stainer       !!
  !!                                               !!
  !!            All rights reserved.               !!
  !!    Licensed under the 3-clause BSD license.   !!
  !!                                               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Module for linked list in fortran
MODULE LINKED_LIST_M
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: INT32
  IMPLICIT NONE
  PRIVATE

  !> custom types
  INTEGER, PARAMETER :: KI4 = INT32

  !> list node type - takes ownership of deallocation of the value pointer in finalize
  TYPE, PUBLIC :: LINKEDLISTNODE
     CLASS(*), POINTER               :: VALUE    => NULL()
     TYPE(LINKEDLISTNODE), POINTER   :: NEXT     => NULL()
     TYPE(LINKEDLISTNODE), POINTER   :: PREV     => NULL()
   CONTAINS
     FINAL :: NODEFINALIZE
  END TYPE LINKEDLISTNODE

  !> list type - takes ownership of deallocation of the value pointer in finalize
  !! Note - do not copy lists, i.e. list1 = list2, this causes memory issues, always pass by reference
  !! Do not return a list from a function
  !! ToDo: Implement copying of lists, do not deep copy pointers
  !! ToDo: Give user option for taking ownership or not (i.e. is user responsible for deallocation of pointers or not)
  TYPE, PUBLIC :: LINKEDLIST
     INTEGER(KI4)                  :: SIZE = 0_KI4
     TYPE(LINKEDLISTNODE), POINTER :: HEAD => NULL()
     TYPE(LINKEDLISTNODE), POINTER :: TAIL => NULL()
   CONTAINS
     PROCEDURE :: APPEND
     PROCEDURE :: REMOVE
     PROCEDURE :: FIRST
     PROCEDURE :: LAST
     PROCEDURE :: ATINDEX
     PROCEDURE :: RESET
     PROCEDURE :: LENGTH
     PROCEDURE :: TRAVERSE
     PROCEDURE, PRIVATE :: CLEANUP
     FINAL :: LISTFINALIZE
  END TYPE LINKEDLIST

  !> extends LinkedList but makes use of a cached last node
  !! if using getatindex in an iterative manner
  !! Improves performance for large lists
  TYPE, EXTENDS(LINKEDLIST), PUBLIC :: CACHEDLINKEDLIST
     PRIVATE
     INTEGER(KI4)                  :: CACHEDLASTINDEX = 0_KI4
     TYPE(LINKEDLISTNODE), POINTER :: CACHEDLASTNODE => NULL()
   CONTAINS
     PROCEDURE :: CACHEDACCESS
     PROCEDURE :: RESET          => RESET_CACHED
     FINAL :: CACHEDLISTFINALIZE
  END TYPE CACHEDLINKEDLIST

CONTAINS

  ! Clean up node - The value is deallocated here
  SUBROUTINE NODEFINALIZE(THIS)
    TYPE(LINKEDLISTNODE), INTENT(INOUT) :: THIS

    IF(ASSOCIATED(THIS%VALUE))THEN
       DEALLOCATE(THIS%VALUE)
       NULLIFY(THIS%VALUE)
       NULLIFY(THIS%NEXT)
       NULLIFY(THIS%PREV)
    END IF
  END SUBROUTINE NODEFINALIZE

  !> Add a value to the list at the tail
  SUBROUTINE APPEND(THIS, VALUE)
    CLASS(LINKEDLIST), INTENT(INOUT) :: THIS
    CLASS(*), INTENT(IN),TARGET      :: VALUE

    TYPE(LINKEDLISTNODE), POINTER :: NODE_PTR, NEXT_PTR, CURRENT_PTR

    ! Create a new node and set the value
    ALLOCATE(NODE_PTR)
    NODE_PTR%VALUE => VALUE
    NODE_PTR%NEXT => NULL()
    THIS%SIZE = THIS%SIZE + 1_KI4

    IF(.NOT. ASSOCIATED(THIS%HEAD))THEN
       THIS%HEAD => NODE_PTR
       THIS%TAIL => NODE_PTR
    ELSE
       THIS%TAIL%NEXT => NODE_PTR
       NODE_PTR%PREV  => THIS%TAIL
       THIS%TAIL      => NODE_PTR
    END IF

  END SUBROUTINE APPEND

  SUBROUTINE REMOVE(THIS, NODE)
    CLASS(LINKEDLIST), INTENT(INOUT)             :: THIS
    TYPE(LINKEDLISTNODE), INTENT(INOUT), POINTER :: NODE

    IF(.NOT. ASSOCIATED(NODE%PREV)) THEN
       THIS%HEAD => NODE%NEXT
       IF(ASSOCIATED(NODE%NEXT)) THEN
          NULLIFY(NODE%NEXT%PREV)
       ENDIF
    ELSE
       NODE%PREV%NEXT => NODE%NEXT
       IF(.NOT. ASSOCIATED(NODE%NEXT)) THEN
          THIS%TAIL => NODE%PREV
       ELSE
          NODE%NEXT%PREV = NODE%PREV
       ENDIF
    ENDIF
    DEALLOCATE(NODE%VALUE)
    NULLIFY(NODE)

  END SUBROUTINE REMOVE

  !> Traverse the list
  SUBROUTINE TRAVERSE(THIS, ITERATOR_FUNC)
    CLASS(LINKEDLIST), INTENT(INOUT) :: THIS
    INTERFACE
       SUBROUTINE ITERATOR_FUNC(NODE)
         IMPORT LINKEDLISTNODE
         TYPE(LINKEDLISTNODE), POINTER, INTENT(INOUT)  :: NODE
       END SUBROUTINE ITERATOR_FUNC
    END INTERFACE

    TYPE(LINKEDLISTNODE), POINTER :: CURRENT_PTR, TEMP_PTR
    INTEGER :: COUNTER

    COUNTER = 0

    CURRENT_PTR => THIS%HEAD
    DO WHILE(ASSOCIATED(CURRENT_PTR) .AND. COUNTER < 10)
       COUNTER = COUNTER + 1
       NULLIFY(TEMP_PTR)
       TEMP_PTR => CURRENT_PTR%NEXT
       CALL ITERATOR_FUNC(CURRENT_PTR)
       CURRENT_PTR => TEMP_PTR
    END DO

  END SUBROUTINE TRAVERSE

  !> Reset the list and cleanup
  SUBROUTINE RESET(THIS)
    CLASS(LINKEDLIST), INTENT(INOUT) :: THIS

    CALL THIS%CLEANUP()

  END SUBROUTINE RESET

  !> Get the size of the list
  PURE FUNCTION LENGTH(THIS) RESULT(SIZE)
    CLASS(LINKEDLIST), INTENT(IN) :: THIS
    INTEGER(KI4) :: SIZE

    SIZE = THIS%SIZE

  END FUNCTION LENGTH

  ! Get the first node
  FUNCTION FIRST(THIS) RESULT(FIRSTNODE)
    CLASS(LINKEDLIST), INTENT(IN) :: THIS
    TYPE(LINKEDLISTNODE), POINTER :: FIRSTNODE

    FIRSTNODE => THIS%HEAD

  END FUNCTION FIRST

  ! Get the last node
  FUNCTION LAST(THIS) RESULT(LASTNODE)
    CLASS(LINKEDLIST), INTENT(IN) :: THIS
    TYPE(LINKEDLISTNODE), POINTER :: LASTNODE

    LASTNODE => THIS%TAIL

  END FUNCTION LAST

  ! Get the node at index
  ! must be between 1 and length()
  FUNCTION ATINDEX(THIS, INDEX) RESULT(INDEXNODE)
    CLASS(LINKEDLIST), INTENT(IN) :: THIS
    INTEGER(KI4), INTENT(IN)      :: INDEX
    TYPE(LINKEDLISTNODE), POINTER :: INDEXNODE

    INTEGER(KI4) :: I

    NULLIFY(INDEXNODE)
    IF(INDEX > 0_KI4 .AND. INDEX <= THIS%SIZE)THEN
       INDEXNODE => THIS%HEAD
       DO I=1, INDEX-1
          INDEXNODE => INDEXNODE%NEXT
       END DO
    END IF

  END FUNCTION ATINDEX

  !> Clean up - deallocation of the nodes in the list
  SUBROUTINE LISTFINALIZE(THIS)
    TYPE(LINKEDLIST), INTENT(INOUT) :: THIS

    CALL THIS%CLEANUP()

  END SUBROUTINE LISTFINALIZE

  !> Clean up - deallocation of the nodes in the list
  SUBROUTINE CLEANUP(THIS)
    CLASS(LINKEDLIST), INTENT(INOUT) :: THIS

    TYPE(LINKEDLISTNODE), POINTER    :: CURRENT_PTR

    CALL THIS%TRAVERSE(DESTROYALL)
    NULLIFY(THIS%HEAD)
    NULLIFY(THIS%TAIL)

  CONTAINS
    SUBROUTINE DESTROYALL(NODE)
      TYPE(LINKEDLISTNODE), POINTER, INTENT(INOUT)  :: NODE

      THIS%HEAD => NODE%NEXT
      DEALLOCATE(NODE)
      NULLIFY(NODE)

      THIS%SIZE = THIS%SIZE - 1_KI4

    END SUBROUTINE DESTROYALL

  END SUBROUTINE CLEANUP

  !> Get the node at index
  ! must be between 1 and length()
  ! It uses the cached index if was set
  ! and then sets the cached node after access
  ! for subsequent calls
  FUNCTION CACHEDACCESS(THIS, INDEX) RESULT(INDEXNODE)
    CLASS(CACHEDLINKEDLIST), INTENT(INOUT) :: THIS
    INTEGER(KI4), INTENT(IN)               :: INDEX

    TYPE(LINKEDLISTNODE), POINTER          :: INDEXNODE

    INTEGER(KI4) :: I
    INTEGER(KI4) :: STARTINDX

    NULLIFY(INDEXNODE)
    IF(INDEX > 0_KI4 .AND. INDEX <= THIS%SIZE)THEN
       ! if last access was cached then use that for speed if we are after it
       IF(THIS%CACHEDLASTINDEX > 0_KI4 .AND. INDEX >= THIS%CACHEDLASTINDEX)THEN
          INDEXNODE => THIS%CACHEDLASTNODE
          STARTINDX =  THIS%CACHEDLASTINDEX
          ! else start at head
       ELSE
          INDEXNODE => THIS%HEAD
          STARTINDX =  1_KI4
       END IF

       DO I=STARTINDX, INDEX-1
          INDEXNODE => INDEXNODE%NEXT
       END DO
       THIS%CACHEDLASTINDEX =  INDEX
       THIS%CACHEDLASTNODE  => INDEXNODE
    END IF

  END FUNCTION CACHEDACCESS

  !> Reset the list and cleanup
  SUBROUTINE RESET_CACHED(THIS)
    CLASS(CACHEDLINKEDLIST), INTENT(INOUT) :: THIS

    CALL THIS%CLEANUP()
    NULLIFY(THIS%CACHEDLASTNODE)
    THIS%CACHEDLASTINDEX = 0_KI4

  END SUBROUTINE RESET_CACHED

  !> Clean up - deallocation of the nodes in the list
  SUBROUTINE CACHEDLISTFINALIZE(THIS)
    TYPE(CACHEDLINKEDLIST), INTENT(INOUT) :: THIS

    CALL THIS%CLEANUP()

  END SUBROUTINE CACHEDLISTFINALIZE

END MODULE LINKED_LIST_M
