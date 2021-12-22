SUBROUTINE MXTRMA(KLX,KVX,KVXS,KIX,PA,PBI,PBS,PX,PY)

!**** *MXTRMA*   - Multiplication of a tridiagonal matrix by a matrix.

!     Purpose.    Multiplication of a tridiagonal matrix by a matrix.
!     --------

!    Example, for KLX=5,KVX=2,KIX=1.

!    I PA (1) PBS(1)   0      0      0    I   I PX(1,1) PX(2,1) I
!    I PBI(1) PA (2) PBS(2)   0      0    I   I PX(1,2) PX(2,2) I
!    I   0    PBI(2) PA (3) PBS(3)   0    I * I PX(1,3) PX(2,3) I
!    I   0      0    PBI(3) PA (4) PBS(4) I   I PX(1,4) PX(2,4) I
!    I   0      0      0    PBI(4) PA (5) I   I PX(1,5) PX(2,5) I

!      I PY(1,1) PY(2,1) I
!      I PY(1,2) PY(2,2) I
!    = I PY(1,3) PY(2,3) I
!      I PY(1,4) PY(2,4) I
!      I PY(1,5) PY(2,5) I

!**   Interface.
!     ----------
!        *CALL* *MXTRMA( ... )

!        Explicit arguments :
!        --------------------
!         KLX:         - Dimension of the matrix.                   (input)
!         KVX,KIX:     - Number of vectors to be multiplied is      (input)
!                        KVX*KIX.
!         KVXS:        - Surdimension corresponding to KVX.         (input)
!         PA:          - Diagonal of the matrix.                    (input)
!         PBI:         - Lower diagonals of the matrix.             (input)
!         PBS:         - Upper diagonals of the matrix.             (input)
!         PX:          - Initial vector:                            (input)
!         PY:          - Final vector:                              (output)

!        Implicit arguments :
!        --------------------
!        none.

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        None.

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        K. YESSAD: FEBRUARY 1995.

!     Modifications.
!     --------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KVXS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KVX 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PA(KLX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBI(KLX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBS(KLX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PX(KVXS,KLX,KIX) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PY(KVXS,KLX,KIX) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JI, JL, JV
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MXTRMA',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1.    COMPUTATION OF PY.
!              ------------------

IF (KVX >= KIX) THEN

  IF (KLX >= 2) THEN

    DO JI=1,KIX
      DO JV=1,KVX
        PY(JV,1,JI) = PA (1)*PX(JV,1,JI)+PBS(1)*PX(JV,2,JI)
      ENDDO
    ENDDO

    DO JI=1,KIX
      DO JV=1,KVX
        DO JL=2,KLX-1
          PY(JV,JL,JI) = PBI(JL-1)*PX(JV,JL-1,JI)&
           & +PA (JL  )*PX(JV,JL  ,JI)&
           & +PBS(JL  )*PX(JV,JL+1,JI)  
        ENDDO
      ENDDO
    ENDDO

    DO JI=1,KIX
      DO JV=1,KVX
        PY(JV,KLX,JI) = PBI(KLX-1)*PX(JV,KLX-1,JI)+PA (KLX  )*PX(JV,KLX  ,JI)
      ENDDO
    ENDDO

  ELSEIF (KLX == 1) THEN

    DO JI=1,KIX
      DO JV=1,KVX
        PY(JV,1,JI) = PA (1)*PX(JV,1,JI)
      ENDDO
    ENDDO

  ENDIF

ELSE

  IF (KLX >= 2) THEN

    DO JV=1,KVX
      DO JI=1,KIX
        PY(JV,1,JI) = PA (1)*PX(JV,1,JI)+PBS(1)*PX(JV,2,JI)
      ENDDO
    ENDDO

    DO JI=1,KIX
      DO JV=1,KVX
        DO JL=2,KLX-1
          PY(JV,JL,JI) = PBI(JL-1)*PX(JV,JL-1,JI)&
           & +PA (JL  )*PX(JV,JL  ,JI)&
           & +PBS(JL  )*PX(JV,JL+1,JI)  
        ENDDO
      ENDDO
    ENDDO

    DO JV=1,KVX
      DO JI=1,KIX
        PY(JV,KLX,JI) = PBI(KLX-1)*PX(JV,KLX-1,JI)+PA (KLX  )*PX(JV,KLX  ,JI)
      ENDDO
    ENDDO

  ELSEIF (KLX == 1) THEN

    DO JV=1,KVX
      DO JI=1,KIX
        PY(JV,1,JI) = PA (1)*PX(JV,1,JI)
      ENDDO
    ENDDO

  ENDIF

ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('MXTRMA',1,ZHOOK_HANDLE)
END SUBROUTINE MXTRMA

