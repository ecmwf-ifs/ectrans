! (C) Copyright 2000- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FOURIER_OUT_MOD
CONTAINS
SUBROUTINE FOURIER_OUT(ZGTF,FOUBUF_IN,KF_FS)

!**** *FOURIER_OUT* - Copy fourier data from local array to buffer

!     Purpose.
!     --------
!        Routine for copying fourier data from local array to buffer

!**   Interface.
!     ----------
!     CALL FOURIER_OUT(...)

!     Explicit arguments :  PREEL - local fourier/GP array
!     --------------------  KFIELDS - number of fields
!
!     Externals.  None.
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 2000-04-01

!     ------------------------------------------------------------------

USE PARKIND_ECTRANS  ,ONLY : JPIM     ,JPRBT

USE TPM_DISTR       ,ONLY : D, MYSETW,D_NPTRLS,D_NSTAGTF,D_NSTAGT0B,D_NPNTGTB0,D_NPROCM,D_NPROCL
USE TPM_GEOMETRY    ,ONLY : G, G_NMEN,G_NLOEN
!

IMPLICIT NONE

REAL(KIND=JPRBT), INTENT(IN) :: ZGTF(:,:)
REAL(KIND=JPRBT), INTENT(OUT) :: FOUBUF_IN(:)
INTEGER(KIND=JPIM),INTENT(IN) :: KF_FS
INTEGER(KIND=JPIM) :: KGL
REAL(KIND=JPRBT)    :: SCAL
INTEGER(KIND=JPIM) :: OFFSET_VAR

INTEGER(KIND=JPIM) :: JM,JF,IGLG,IPROC,IR,II,ISTA, ISTA1,JMMAX, iunit

INTEGER(KIND=JPIM) :: IBEG,IEND,IINC, IOFF,iimax1,iimax2,iimax3

! scale results and move into next transformation buffer

OFFSET_VAR=D_NPTRLS(MYSETW)
!$ACC DATA PRESENT(D,G_NLOEN,D_NSTAGT0B,D_NPNTGTB0,FOUBUF_IN,D_NPROCM,ZGTF,G_NMEN,D_NSTAGTF)
!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(IGLG,IOFF,SCAL,JM) DEFAULT(NONE)
DO KGL=1,D%NDGL_FS
  DO JF=1,KF_FS
    IGLG = OFFSET_VAR+KGL-1
    IOFF = D_NSTAGTF(KGL)+1

    SCAL = 1._JPRBT/REAL(G_NLOEN(IGLG),JPRBT)
    !$ACC LOOP SEQ
    DO JM=0,G_NMEN(IGLG)
      IPROC = D_NPROCM(JM)
      ISTA  = (D_NSTAGT0B(IPROC)+D_NPNTGTB0(JM,KGL))*2*KF_FS

      FOUBUF_IN(ISTA+2*JF-1) = SCAL * ZGTF(2*JF-1, 2*JM+IOFF)
      FOUBUF_IN(ISTA+2*JF  ) = SCAL * ZGTF(2*JF  , 2*JM+IOFF)
     ENDDO
   ENDDO
ENDDO
!$ACC END DATA

!     ------------------------------------------------------------------

END SUBROUTINE FOURIER_OUT
END MODULE FOURIER_OUT_MOD

