! (C) Copyright 2000- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE VDTUV_MOD
CONTAINS
SUBROUTINE VDTUV(KFIELD,PEPSNM,PVOR,PDIV,PU,PV)

USE PARKIND_ECTRANS  ,ONLY : JPIM     ,JPRB,  JPRBT

USE TPM_DIM         ,ONLY : R, R_NTMAX
USE TPM_FIELDS      ,ONLY : F
USE TPM_DISTR       ,ONLY : D,D_NUMP,D_MYMS
use tpm_gen, only: nout


!**** *VDTUV* - Compute U,V in  spectral space

!     Purpose.
!     --------
!        In Laplace space compute the the winds
!        from vorticity and divergence.

!**   Interface.
!     ----------
!        CALL VDTUV(...)

!        Explicit arguments :  KM -zonal wavenumber (input-c)
!        --------------------  KFIELD - number of fields (input-c)
!                              PEPSNM - REPSNM for wavenumber KM (input-c)
!                              PVOR(NLEI1,2*KFIELD) - vorticity (input)
!                              PDIV(NLEI1,2*KFIELD) - divergence (input)
!                              PU(NLEI1,2*KFIELD)   - u wind (output)
!                              PV(NLEI1,2*KFIELD)   - v wind (output)
!        Organisation within NLEI1:
!        NLEI1 = NSMAX+4+mod(NSMAX+4+1,2)
!                        overdimensioning
!        1        : n=NSMAX+2
!        2        : n=NSMAX+1
!        3        : n=NSMAX
!        .        :
!        .        :
!        NSMAX+3  : n=0
!        NSMAX+4  : n=-1

!        Implicit arguments :  Eigenvalues of inverse Laplace operator
!        --------------------  from YOMLAP

!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Temperton, 1991, MWR 119 p1303

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From VDTUV in IFS CY22R1

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM) :: KM, kmloc
INTEGER(KIND=JPIM), INTENT(IN) :: KFIELD
REAL(KIND=JPRBT), INTENT(IN)    :: PEPSNM(1:D%NUMP,0:R%NTMAX+2)
REAL(KIND=JPRB), INTENT(INOUT)    :: PVOR(:,:,:),PDIV(:,:,:)
REAL(KIND=JPRB), INTENT(OUT)   :: PU  (:,:,:),PV  (:,:,:)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: II, IJ, IR, J, JN, JI

!     LOCAL REAL SCALARS
REAL(KIND=JPRBT) :: ZKM

!$ACC DATA                                     &
!$ACC      COPYIN (D,D%MYMS,F,F%RLAPIN,F%RN)   &
!$ACC      PRESENT(PEPSNM, PVOR, PDIV)         &
!$ACC      PRESENT(PU, PV, D_MYMS)


!     ------------------------------------------------------------------

!*       1.    COMPUTE U V FROM VORTICITY AND DIVERGENCE.
!              ------------------------------------------

!$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE)
DO KMLOC=1,D%NUMP
  DO J=1,KFIELD
    IR = 2*J-1
    II = IR+1
    KM = D_MYMS(KMLOC)
    ZKM = REAL(KM,JPRBT)

    IF(KM == 0) THEN
      !$ACC LOOP SEQ
      DO JN=0,R_NTMAX+1
        JI = R_NTMAX+3-JN
        PU(IR,JI,KMLOC) = +&
         &(JN-1)*PEPSNM(KMLOC,JN)*F%RLAPIN(JN-1)*PVOR(IR,JI+1,KMLOC)-&
         &(JN+2)*PEPSNM(KMLOC,JN+1)*F%RLAPIN(JN+1)*PVOR(IR,JI-1,KMLOC)
        PV(IR,JI,KMLOC) = -&
         &(JN-1)*PEPSNM(KMLOC,JN)*F%RLAPIN(JN-1)*PDIV(IR,JI+1,KMLOC)+&
         &(JN+2)*PEPSNM(KMLOC,JN+1)*F%RLAPIN(JN+1)*PDIV(IR,JI-1,KMLOC)
      ENDDO

    ELSE
      !$ACC LOOP SEQ
      DO JN=KM,R_NTMAX+1
        JI = R_NTMAX+3-JN
        PU(ir,JI,kmloc) = -ZKM*F%RLAPIN(JN)*PDIV(ii,JI,kmloc)+&
         &(JN-1)*PEPSNM(KMLOC,JN)*F%RLAPIN(JN-1)*PVOR(ir,JI+1,kmloc)-&
         &(JN+2)*PEPSNM(KMLOC,JN+1)*F%RLAPIN(JN+1)*PVOR(ir,JI-1,kmloc)
        PU(ii,JI,kmloc) = +ZKM*F%RLAPIN(JN)*PDIV(ir,JI,kmloc)+&
         &(JN-1)*PEPSNM(KMLOC,JN)*F%RLAPIN(JN-1)*PVOR(ii,JI+1,kmloc)-&
         &(JN+2)*PEPSNM(KMLOC,JN+1)*F%RLAPIN(JN+1)*PVOR(ii,JI-1,kmloc)
        PV(ir,JI,kmloc) = -ZKM*F%RLAPIN(JN)*PVOR(ii,JI,kmloc)-&
         &(JN-1)*PEPSNM(KMLOC,JN)*F%RLAPIN(JN-1)*PDIV(ir,JI+1,kmloc)+&
         &(JN+2)*PEPSNM(KMLOC,JN+1)*F%RLAPIN(JN+1)*PDIV(ir,JI-1,kmloc)
        PV(ii,JI,kmloc) = +ZKM*F%RLAPIN(JN)*PVOR(ir,JI,kmloc)-&
         &(JN-1)*PEPSNM(KMLOC,JN)*F%RLAPIN(JN-1)*PDIV(ii,JI+1,kmloc)+&
         &(JN+2)*PEPSNM(KMLOC,JN+1)*F%RLAPIN(JN+1)*PDIV(ii,JI-1,kmloc)
      ENDDO
    ENDIF
  ENDDO
ENDDO

!$ACC END DATA
!     ------------------------------------------------------------------

END SUBROUTINE VDTUV
END MODULE VDTUV_MOD

