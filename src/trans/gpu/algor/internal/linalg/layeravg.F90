!*********************** SUBROUTINE LAYERAVG *************************

   SUBROUTINE LAYERAVG(LDGRADPS,PX1,PX2,PY2,KN1,KN2,KI,PZ,PZS,PZPS)

!-----------------------------------------------------------------------

!     This software was developed and tested by Y.J. Rochon
!     with contributions from L. Garand, D.S. Turner, and S. Polavarapu
!     of the Atmospheric Science and Technology Directorate,
!     Environment Canada (EC).

!     Copyright 2006, Environment Canada, All Rights Reserved.

!     Please send any improvements to any one of the above authors
!     at firstname.lastname@ec.gc.ca (e.g yves.rochon@ec.gc.ca) and provide
!     appropriate aknowledgement through referencing to the authors or
!     the corresponding journal article (see below).

!     Description:

!     Perform integration between points PX1(KI-1) and
!     PX1(KI+1) with  piecewise linear weighting having
!     weights of zero at  ki-1 and ki+1 and max weight at ki.

!     Method:
!     - Piecewise weighted interpolation in ln(P).
!     - Journal reference:
!       Rochon, Y.J., L. Garand, D.S. Turner, and S. Polavarapu.
!       Jacobian mapping between coordinate systems in data assimilation,
!       Q. J. of the Royal Met. Soc., Accepted pending revisions on 25
!       Sept. 2006.

!     History:

!     Version   Date        Comment
!     -------   ----        -------
!     1         10/2005     Original F77 code by Y.J. Rochon with
!                           contributions from S. Polavarapu. (EC)
!     2         10/2006     Adaptation for F90 and for improved consistency
!                           with RTTOV code format (Y.J. Rochon, EC).
!                           Argument names of original code were not changed.
!     
!     Arguments:
!     
!     Input

!     PX1(KN1).................Reference levels (e.g. lnP; in increasing values)
!     PX2(KN2).................Available levels (e.g. lnP; in increasing values)
!     PY2(KN2).................Values at PX2 levels (e.g. temperature)
!     LDGRADPS..................Flag indicating if gradient w.r.t. Ps is 
!                              required. True for yes.
!     PZS......................dlnP/dPs required when  LDGRADPS=.true.

!     KN1......................Dimension of PX1.
!     KN2......................Dimension of other arrays.
!     KI.......................Identifies region of interest: 
!                              PX1(KI-1) to PX1(KI+1)   
!     LDGRADPS..................Flag indicating of gradient w.r.t. Ps is 
!                              required. True for yes.
!     PZ.......................Extra array related to gradients w.r.t. Ps
!                              when LDGRADPS=.true.
!     PZS......................dlnP/dPs required when  LDGRADPS=.true.

!     Output

!     PZ(KN2)..................Resultant accumulated contribution factors
!                              (a row of the K matrix)
!     PZPS.....................Gradient term w.r.t. Ps when LDGRADPS=.true.
!                              (extra element to the row of the K matrix)

!     External functions and subroutines:
!     
!     LAYERAVG: Performs integration between points ZLNPO(JO-1) and
!               ZLNPO(JO+1) with  piecewise linear weighting having
!               weights of zero at  JO-1 and JO+1 and max weight at JO.

!     Related routines: INTAVGAD, INTAVG, INTAVGTL, INTAVGK, SUBLAYER

!     Assumptions:

!     1) PX1(i)<PX1(i+1) & PX2(i)<PX2(i+1)

!     Comments:
!                    
!     1) The gradient calc w.r.t Ps does not account for possible
!     re-assigment of lowest contributing RTM level in INTAVG* and LAYERAVG
!     to the surface model level when the original values are below
!     the model surface. This will result in differences from gradients
!     calculated using finite differences.

!--------------------------------------------------------------------
!     
      USE PARKIND1  ,ONLY : JPIM     ,JPRB
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

      IMPLICIT NONE
!     
! --- Subroutine arguments

      LOGICAL, INTENT(IN) :: LDGRADPS                  
      INTEGER(KIND=JPIM), INTENT(IN) :: KN1,KN2,KI
      REAL(KIND=JPRB), INTENT(IN)    :: PX1(KN1),PX2(KN2),PY2(KN2),PZS(KN2)
      REAL(KIND=JPRB), INTENT(INOUT) :: PZ(KN2),PZPS
      REAL(KIND=JPRB)    :: ZHOOK_HANDLE
!                  
! --- Local scalars and arrays
!                          
      INTEGER(KIND=JPIM) :: J,IC,ISKIP
      REAL(KIND=JPRB)    :: Z1,Z2,Z3,ZW1,ZW2,ZSUM


#include "sublayer.h"
!-----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('LAYERAVG',0,ZHOOK_HANDLE)

! --- Identify boundary points

      Z2=PX1(KI)

      IF (KI == 1) THEN 
         Z1=2.0*Z2-PX1(KI+1)
      ELSE
         Z1=PX1(KI-1)
      ENDIF   

      IF (KI == KN1) THEN
         Z3=2.0*Z2-Z1
      ELSE   
         Z3=PX1(KI+1)
      ENDIF
      IF (Z3 > PX2(KN2)) Z3=PX2(KN2)

      ISKIP=0
      IF (Z2 >= PX2(KN2)) THEN
         Z3=PX2(KN2)
         Z2=PX2(KN2)
         ISKIP=1
      ENDIF

! --- Determine forward interpolator (kflag=0) or TLM (kflag>0)

      PZPS=0.0
      PZ(1:KN2)=0.0
      IC=0
      DO J=1,KN2-1
         IF (PX2(J) >= Z3) GO TO 1000

         IF (PX2(J) <= Z2.AND.PX2(J+1) > Z1) THEN 

            CALL SUBLAYER(Z1,Z2,Z3,PX2(J),PX2(J+1),LDGRADPS, &
                      &    PY2(J),PY2(J+1),ZW1,ZW2,          &
                      &    PZS(J),PZS(J+1),PZPS)
            PZ(J)=PZ(J)+ZW1
            PZ(J+1)=PZ(J+1)+ZW2
            IC=1
         ENDIF

         IF (PX2(J) < Z3.AND.PX2(J+1) >= Z2.AND.ISKIP == 0) THEN

            CALL SUBLAYER(Z3,Z2,Z1,PX2(J),PX2(J+1),LDGRADPS, &
                      &    PY2(J),PY2(J+1),ZW1,ZW2,          &
                      &    PZS(J),PZS(J+1),PZPS)
            PZ(J)=PZ(J)+ZW1
            PZ(J+1)=PZ(J+1)+ZW2
            IC=1
         ENDIF
      ENDDO
      J=KN2
 1000 CONTINUE
      IF (IC == 0) PZ(J)=1.0

!     Normalize sum to unity (instead of calculating and dividing by
!     weighting denominator)

      ZSUM=SUM(PZ(1:KN2))
      PZ(1:KN2)=PZ(1:KN2)/ZSUM
      IF (IC /= 0) PZPS=PZPS/ZSUM

   IF (LHOOK) CALL DR_HOOK('LAYERAVG',1,ZHOOK_HANDLE)
   END SUBROUTINE LAYERAVG
