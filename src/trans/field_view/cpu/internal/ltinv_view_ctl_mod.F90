! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LTINV_VIEW_CTL_MOD
CONTAINS
SUBROUTINE LTINV_VIEW_CTL(KF_OUT_LT,&
 & YDSPVVOR,YDSPVDIV,YDSPVSCALAR,&
 & FSPGL_PROC)

!**** *LTINV_VIEW_CTL* - Control routine for inverse Legandre transform.

!     Purpose.
!     --------
!        Control routine for the inverse LEGENDRE transform

!**   Interface.
!     ----------
!     CALL INV_TRANS_CTL(...)
!     KF_OUT_LT    - number of fields coming out from inverse LT
!     PSPVOR(:,:)  - spectral vorticity (input)
!     PSPDIV(:,:)  - spectral divergence (input)
!     PSPSCALAR(:,:) - spectral scalarvalued fields (input)
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

USE PARKIND1   ,ONLY : JPIM     ,JPRB

USE TPM_GEN    ,ONLY : LALLOPERM
USE TPM_TRANS  ,ONLY : FOUBUF, FOUBUF_IN, LSCDERS
USE TPM_DISTR  ,ONLY : D
USE TPM_FLT    ,ONLY : S

USE LTINV_VIEW_MOD,ONLY : LTINV_VIEW
USE TRMTOL_MOD ,ONLY : TRMTOL
USE ECTRANS_FIELD_VIEW_INTERNAL_UTIL_MOD, ONLY: SPEC_VIEW

IMPLICIT NONE

#include "fspgl_intf.h"

INTEGER(KIND=JPIM),INTENT(IN) :: KF_OUT_LT
TYPE(SPEC_VIEW), INTENT(IN) :: YDSPVVOR(:), YDSPVDIV(:)
TYPE(SPEC_VIEW), INTENT(IN) :: YDSPVSCALAR(:)

PROCEDURE (FSPGL_INTF), POINTER, OPTIONAL, INTENT(IN)  :: FSPGL_PROC

INTEGER(KIND=JPIM) :: JM,IM,IBLEN,ILEI2,IDIM1
INTEGER(KIND=JPIM) :: IF_UV,IF_SCALARS,IF_SCDERS
!     IF_UV        - local number of spectral u-v fields
!     IF_SCALARS   - local number of scalar spectral fields
!     IF_SCDERS    - local number of derivatives of scalar spectral fields

!     ------------------------------------------------------------------

CALL GSTATS(102,0)

IF_UV = SIZE(YDSPVVOR)
IF_SCALARS = SIZE(YDSPVSCALAR)
IF_SCDERS = 0
IF (LSCDERS)  IF_SCDERS = SIZE(YDSPVSCALAR)

ILEI2 = 8*IF_UV + 2*IF_SCALARS + 2*IF_SCDERS
IDIM1 = 2*KF_OUT_LT
IBLEN = D%NLENGT0B*2*KF_OUT_LT
IF (ALLOCATED(FOUBUF)) THEN
  IF (MAX(1,IBLEN) > SIZE(FOUBUF)) THEN
    DEALLOCATE(FOUBUF)
    ALLOCATE(FOUBUF(MAX(1,IBLEN)))
  ENDIF
ELSE
  ALLOCATE(FOUBUF(MAX(1,IBLEN)))
ENDIF
IF (ALLOCATED(FOUBUF_IN)) THEN
  IF (MAX(1,IBLEN) > SIZE(FOUBUF_IN)) THEN
    DEALLOCATE(FOUBUF_IN)
    ALLOCATE(FOUBUF_IN(MAX(1,IBLEN)))
  ENDIF
ELSE
  ALLOCATE(FOUBUF_IN(MAX(1,IBLEN)))
  FOUBUF_IN(:) = 0
ENDIF

! Following switch necessary when latlon grids are used with different increments in NS and EW direction.
! Otherwise unassigned values will appear in output. This is very likely a bug (ATLAS-149)
IF (S%LDLL) THEN
  FOUBUF_IN(:) = 0
ENDIF

IF(KF_OUT_LT > 0) THEN
  CALL GSTATS(1647,0)

  !!!WARNING!!! Duplication of code besides the FSPGL_PROC argument.
              ! It seems that gfortran 10 does not retain the value
              ! of FSPGL_PROC within the OMP region.
  IF( PRESENT(FSPGL_PROC) ) THEN
    !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JM,IM)
    DO JM=1,D%NUMP
      IM = D%MYMS(JM)
      CALL LTINV_VIEW(IM,JM,KF_OUT_LT,IF_UV,IF_SCALARS,IF_SCDERS,ILEI2,IDIM1,&
       & YDSPVVOR,YDSPVDIV,YDSPVSCALAR,&
       &FSPGL_PROC)
    ENDDO
    !$OMP END PARALLEL DO
  ELSE
    !$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JM,IM)
    DO JM=1,D%NUMP
      IM = D%MYMS(JM)
      CALL LTINV_VIEW(IM,JM,KF_OUT_LT,IF_UV,IF_SCALARS,IF_SCDERS,ILEI2,IDIM1,&
       & YDSPVVOR,YDSPVDIV,YDSPVSCALAR)
    ENDDO
    !$OMP END PARALLEL DO
  ENDIF
  CALL GSTATS(1647,1)
ENDIF

CALL GSTATS(102,1)

CALL GSTATS(152,0)
CALL TRMTOL(FOUBUF_IN,FOUBUF,2*KF_OUT_LT)
CALL GSTATS(152,1)
IF (.NOT.LALLOPERM) DEALLOCATE(FOUBUF_IN)
!     ------------------------------------------------------------------

END SUBROUTINE LTINV_VIEW_CTL
END MODULE LTINV_VIEW_CTL_MOD
