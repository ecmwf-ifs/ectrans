MODULE LTINV_CTL_MOD
CONTAINS
SUBROUTINE LTINV_CTL(KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS,&
 & PSPVOR,PSPDIV,PSPSCALAR,KFLDPTRUV,KFLDPTRSC,FSPGL_PROC)

!**** *LTINV_CTL* - Control routine for inverse Legandre transform.

!     Purpose.
!     --------
!        Control routine for the inverse LEGENDRE transform

!**   Interface.
!     ----------
!     CALL INV_TRANS_CTL(...)
!     KF_OUT_LT    - number of fields coming out from inverse LT
!     KF_UV        - local number of spectral u-v fields
!     KF_SCALARS   - local number of scalar spectral fields
!     KF_SCDERS    - local number of derivatives of scalar spectral fields
!     PSPVOR(:,:)  - spectral vorticity (input)
!     PSPDIV(:,:)  - spectral divergence (input)
!     PSPSCALAR(:,:) - spectral scalarvalued fields (input)
!     KFLDPTRUV(:) - field pointer array for vor./div.
!     KFLDPTRSC(:) - field pointer array for PSPSCALAR
!     FSPGL_PROC  - external procedure to be executed in fourier space
!                   before transposition  

!     Method.
!     -------

!     Externals.  
!     ----------  
!                 

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-06-03

!     ------------------------------------------------------------------

#include "tsmbkind.h"

USE TPM_GEN
USE TPM_DIM
USE TPM_TRANS
USE TPM_DISTR

USE LTINV_MOD
USE TRMTOL_MOD

IMPLICIT NONE

INTEGER_M,INTENT(IN) :: KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS
REAL_B ,OPTIONAL, INTENT(IN)  :: PSPVOR(:,:)
REAL_B ,OPTIONAL, INTENT(IN)  :: PSPDIV(:,:)
REAL_B ,OPTIONAL, INTENT(IN)  :: PSPSCALAR(:,:)
INTEGER_M,OPTIONAL,INTENT(IN) :: KFLDPTRUV(:)
INTEGER_M,OPTIONAL,INTENT(IN) :: KFLDPTRSC(:)
EXTERNAL  FSPGL_PROC
OPTIONAL  FSPGL_PROC

INTEGER_M :: JM,IM,IBLEN,ILEI2,IDIM1

!     ------------------------------------------------------------------

ILEI2 = 8*KF_UV + 2*KF_SCALARS + 2*KF_SCDERS
IDIM1 = 2*KF_OUT_LT
IBLEN = D%NLENGT0B*2*KF_OUT_LT
IF(LALLOPERM) THEN
  IF(ALLOCATED(FOUBUF_IN)) THEN
    IF(SIZE(FOUBUF_IN) < IBLEN) THEN
      DEALLOCATE(FOUBUF_IN)
    ENDIF
  ENDIF
ENDIF
IF(.NOT. ALLOCATED(FOUBUF_IN)) ALLOCATE(FOUBUF_IN(IBLEN))

CALL GSTATS(102,0)
IF(KF_OUT_LT > 0) THEN
#ifndef HLOMP
!$OMP PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JM,IM)
#endif
  DO JM=1,D%NUMP
    IM = D%MYMS(JM)
    CALL LTINV(IM,JM,KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS,ILEI2,IDIM1,&
     & PSPVOR,PSPDIV,PSPSCALAR,KFLDPTRUV,KFLDPTRSC,FSPGL_PROC)
  ENDDO
#ifndef HLOMP
!$OMP END PARALLEL DO
#endif
ENDIF
CALL GSTATS(102,1)

IBLEN = D%NLENGT0B*2*KF_OUT_LT
IF(LALLOPERM) THEN
  IF(ALLOCATED(FOUBUF)) THEN
    IF(SIZE(FOUBUF) < IBLEN) THEN
      DEALLOCATE(FOUBUF)
    ENDIF
  ENDIF
ENDIF
IF(.NOT. ALLOCATED(FOUBUF))  ALLOCATE(FOUBUF(IBLEN))

CALL GSTATS(152,0)
CALL TRMTOL(FOUBUF_IN,FOUBUF,2*KF_OUT_LT)
CALL GSTATS(152,1)
IF(.NOT. LALLOPERM) THEN
  DEALLOCATE(FOUBUF_IN)
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE LTINV_CTL
END MODULE LTINV_CTL_MOD
