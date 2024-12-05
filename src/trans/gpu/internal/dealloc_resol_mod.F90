! (C) Copyright 2013- ECMWF.
! (C) Copyright 2013- Meteo-France.
! (C) Copyright 2024- NVIDIA.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE DEALLOC_RESOL_MOD
CONTAINS
SUBROUTINE DEALLOC_RESOL(KRESOL)

!**** *DEALLOC_RESOL* - Deallocations of a resolution

!     Purpose.
!     --------
!     Release allocated arrays for a given resolution

!**   Interface.
!     ----------
!     CALL DEALLOC_RESOL

!     Explicit arguments : KRESOL : resolution tag
!     --------------------

!     Method.
!     -------

!     Externals.  None
!     ----------

!     Author.
!     -------
!        R. El Khatib *METEO-FRANCE*

!     Modifications.
!     --------------
!        Original : 09-Jul-2013 from trans_end

!     ------------------------------------------------------------------

USE PARKIND_ECTRANS, ONLY: JPIM
USE TPM_DIM,         ONLY: R, DIM_TYPE
USE TPM_GEN,         ONLY: LENABLED, NOUT, NDEF_RESOL
USE TPM_DISTR,       ONLY: D, DISTR_TYPE, NPRTRV
USE TPM_GEOMETRY,    ONLY: G, GEOM_TYPE
USE TPM_FIELDS,      ONLY: F, FIELDS_TYPE
USE TPM_FIELDS_GPU,  ONLY: FG, FIELDS_GPU_TYPE
USE TPM_HICFFT,      ONLY: CLEAN_FFT
USE HICBLAS_MOD,     ONLY: CLEAN_GEMM
USE TPM_FLT,         ONLY: S, FLT_TYPE_WRAP
USE TPM_CTL,         ONLY: C
USE SEEFMM_MIX,      ONLY: FREE_SEEFMM
USE SET_RESOL_MOD,   ONLY: SET_RESOL
!

IMPLICIT NONE

INTEGER(KIND=JPIM),  INTENT(IN) :: KRESOL
INTEGER(KIND=JPIM) :: JMLOC,IPRTRV,JSETV,IMLOC,IM,ILA,ILS, JRESOL
TYPE(DIM_TYPE) :: R_
TYPE(DISTR_TYPE) :: D_
TYPE(GEOM_TYPE) :: G_
TYPE(FIELDS_TYPE) :: F_
TYPE(FIELDS_GPU_TYPE) :: FG_
TYPE(FLT_TYPE_WRAP) :: S_

!     ------------------------------------------------------------------

IF (.NOT.LENABLED(KRESOL)) THEN

  WRITE(UNIT=NOUT,FMT='('' DEALLOC_RESOL WARNING : KRESOL = '',I3,'' ALREADY DISABLED '')') KRESOL

ELSE

  CALL SET_RESOL(KRESOL)

#ifdef ACCGPU
!$ACC EXIT DATA DELETE(R) ASYNC(1)
!$ACC EXIT DATA DELETE(FG,FG%ZAA0,FG%ZAS0) IF(ALLOCATED(FG%ZAA0)) ASYNC(1)
!$ACC EXIT DATA DELETE(FG,FG%ZAA,FG%ZAS,FG%ZEPSNM) ASYNC(1)
!$ACC EXIT DATA DELETE(F,F%RLAPIN,F%RACTHE,F%RW) ASYNC(1)
!$ACC EXIT DATA DELETE(D,D%MYMS,D%NPNTGTB0,D%NPNTGTB1,D%NSTAGT0B,D%NSTAGT1B,D%NSTAGTF,D%NPROCM)&
!$ACC&          DELETE(D%NPTRLS,D%MSTABF,D%NASM0,D%OFFSETS_GEMM1,D%OFFSETS_GEMM2) ASYNC(1)
!$ACC EXIT DATA DELETE(G,G%NDGLU,G%NMEN,G%NLOEN) ASYNC(1)
!$ACC WAIT(1)
#endif
#ifdef OMPGPU
#endif

  ! TPM_FLD is more complex because it has pointers
  IF( ALLOCATED(S%FA) ) THEN
    DO JMLOC=1,D%NUMP,NPRTRV  ! +++++++++++++++++++++ JMLOC LOOP ++++++++++
      IPRTRV=MIN(NPRTRV,D%NUMP-JMLOC+1)
      DO JSETV=1,IPRTRV
        IMLOC=JMLOC+JSETV-1
        IM = D%MYMS(IMLOC)
        ILA = (R%NSMAX-IM+2)/2
        ILS = (R%NSMAX-IM+3)/2
        IF(.NOT.C%CIO_TYPE == 'mbuf' .AND. ASSOCIATED(S%FA(IMLOC)%RPNMA)) DEALLOCATE(S%FA(IMLOC)%RPNMA)
        IF(.NOT.C%CIO_TYPE == 'mbuf' .AND. ASSOCIATED(S%FA(IMLOC)%RPNMS)) DEALLOCATE(S%FA(IMLOC)%RPNMS)
        IF(S%LDLL) THEN
          IF(.NOT.C%CIO_TYPE == 'mbuf' .AND. ASSOCIATED(S%FA(IMLOC)%RPNMWI)) DEALLOCATE(S%FA(IMLOC)%RPNMWI)
          IF(.NOT.C%CIO_TYPE == 'mbuf' .AND. ASSOCIATED(S%FA(IMLOC)%RPNMWO)) DEALLOCATE(S%FA(IMLOC)%RPNMWO)
        ENDIF
      ENDDO
    ENDDO
    DEALLOCATE(S%FA)
  ENDIF
  IF(S%LDLL) THEN
    CALL FREE_SEEFMM(S%FMM_INTI)
    IF(ASSOCIATED(S%FMM_INTI)) DEALLOCATE(S%FMM_INTI)
  ENDIF
  S = S_

  ! Empty all fields (none of them has pointers; allocatable arrays implicitly deallocate)
  D = D_
  F = F_
  FG = FG_
  R = R_
  G = G_

  CALL CLEAN_FFT(KRESOL)
  CALL CLEAN_GEMM(KRESOL)

  LENABLED(KRESOL)=.FALSE.
  NDEF_RESOL = COUNT(LENABLED)
  ! Do not stay on a disabled resolution
  DO JRESOL=1,SIZE(LENABLED)
    IF (LENABLED(JRESOL)) THEN
      CALL SET_RESOL(JRESOL)
      EXIT
    ENDIF
  ENDDO

ENDIF
!     ------------------------------------------------------------------

END SUBROUTINE DEALLOC_RESOL
END MODULE DEALLOC_RESOL_MOD
