SUBROUTINE MXMAOP(PA,KA,KAD,PB,KB,KBD,PC,KC,KCA,KAR,KAC,KBC)

!**** *MXMAOP - Optimize call to SGEMMX

!     Purpose.
!     --------
!        Make sure SGEMMX is called in a way to insure maximum optimization.

!**   Interface.
!     ----------
!        CALL MXMAOP(PA,KA,KAD,PB,KB,KBD,PC,KC,KCA,KAR,KAC,KBC)

!        Explicit arguments : See SGEMMX documentaion.
!        --------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals.   SGEMMX in Cray library.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 88-01-28
!     ------------------------------------------------------------------

#include "tsmbkind.h"

IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KA
INTEGER_M :: KAC
INTEGER_M :: KAD
INTEGER_M :: KAR
INTEGER_M :: KB
INTEGER_M :: KBC
INTEGER_M :: KBD
INTEGER_M :: KC
INTEGER_M :: KCA


REAL_B :: PA(1),PB(1),PC(2)

!     ------------------------------------------------------------------

!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------
IF (KAR >= KBC) THEN
  CALL SGEMMX(KAR,KBC,KAC,_ONE_,PA,KA,KAD,PB,KB,KBD,_ZERO_,PC,KC,KCA)
ELSE
  CALL SGEMMX(KBC,KAR,KAC,_ONE_,PB,KBD,KB,PA,KAD,KA,_ZERO_,PC,KCA,KC)
ENDIF

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE MXMAOP



