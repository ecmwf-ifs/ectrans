MODULE PE2SET_MOD
CONTAINS
SUBROUTINE PE2SET(KPE,KPRGPNS,KPRGPEW,KPRTRW,KPRTRV)

#ifdef DOC

!**** *PE2SET* - Convert from PE number to set numbers

!     Purpose.
!     --------
!        Convert from PE number to set numbers in both
!                  grid-point space and spectral space

!**   Interface.
!     ----------
!        *CALL* *PE2SET(KPE,KPRGPNS,KPRGPEW,KPRTRW,KPRTRV)

!        Explicit arguments :  
!        --------------------
!                  input:   KPE     - integer processor number 
!                                     in the range 1 .. NPROC
!                  output:  KPRGPNS - integer A set number in grid space
!                                     in the range 1 .. NPRGPNS
!                           KPRGPEW - integer B set number in grid space
!                                     in the range 1 .. NPRGPEW
!                           KPRTRW  - integer A set number in spectral space
!                                     in the range 1 .. NPRTRW 
!                           KPRTRV  - integer B set number in spectral space
!                                     in the range 1 .. NPRTRV 

!        Implicit arguments :  YOMMP parameters
!                              NPRGPNS,NPRGPEW,NPRTRW,NPRTRV,NPROC

!        --------------------
!     Method.
!     -------

!        PE allocation order is row oriented (e.g. NPRGPNS or NPRTRW = 4):

!                1  2  3  4 
!                5  6  7  8 
!                9 10 11 12 
!               13 14 15 16 
!                .  .  .  .

!     Externals.
!     ----------
!         NONE

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        David Dent *ECMWF*

!     Modifications.
!     --------------
!        Original : 98-08-19
!        Revision : 98-10-13 row ordering
!     ------------------------------------------------------------------
#endif

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DISTR
USE ABORT_TRANS_MOD


IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)  :: KPE
INTEGER(KIND=JPIM),INTENT(OUT) :: KPRGPNS,KPRGPEW,KPRTRW,KPRTRV

!     ------------------------------------------------------------------

!*       1.    Check input argument for validity 
!              ---------------------------------

IF(KPE <= 0.OR.KPE > NPROC) THEN
  WRITE(*,'(A,2I8)') ' PE2SET INVALID ARGUMENT ',KPE,NPROC
  CALL ABORT_TRANS(' PE2SET INVALID ARGUMENT ')

ELSE

!*       2.    Compute output parameters
!              -------------------------

  KPRGPEW=MOD(KPE-1,NPRGPEW)+1
  KPRGPNS=(KPE-1)/NPRGPEW+1
  KPRTRV =MOD(KPE-1,NPRTRV)+1
  KPRTRW =(KPE-1)/NPRTRV+1

ENDIF

END SUBROUTINE PE2SET
END MODULE PE2SET_MOD
