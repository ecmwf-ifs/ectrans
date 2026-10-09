! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE INV_TRANS_VIEW_CTL_MOD
CONTAINS
  SUBROUTINE INV_TRANS_VIEW_CTL(KPROMA,KGPBLKS, &
                            & YDSPVVOR, YDSPVDIV, YDSPVSCALAR, &
                            & YDGVU,YDGVV,&
                            & YDGVVOR,YDGVDIV,&
                            & YDGVSCALAR,&
                            & YDGVU_EW,YDGVV_EW,&
                            & YDGVSCALAR_NS, YDGVSCALAR_EW,&
                            & FSPGL_PROC)

    !**** *INV_TRANS_VIEW_CTL* - Control routine for inverse spectral transform.

    !     Purpose.
    !     --------
    !        Control routine for the inverse spectral transform

    !**   Interface.
    !     ----------
    !     CALL INV_TRANS_VIEW_CTL(...)

    !     Explicit arguments :
    !     --------------------
    
    !     KF_OUT_LT    - total number of fields coming out from inverse LT
    !     IF_UV        - local number of spectral u-v fields
    !     PSPVOR(:,:)  - spectral vorticity (input)
    !     PSPDIV(:,:)  - spectral divergence (input)
    !     PSPSCALAR(:,:) - spectral scalarvalued fields (input)
    !     KVSETUV(:)  - indicating which 'b-set' in spectral space owns a
    !                   vor/div field. Equivalant to NBSETLEV in the IFS.
    !                   The length of KVSETUV should be the GLOBAL number
    !                   of u/v fields which is the dimension of u and v releated
    !                   fields in grid-point space.
    !     KVESETSC(:) - indicating which 'b-set' in spectral space owns a
    !                   scalar field. As for KVSETUV this argument is required
    !                   if the total number of processors is greater than
    !                   the number of processors used for distribution in
    !                   spectral wave space.
    !     FSPGL_PROC  - external procedure to be executed in fourier space
    !                   before transposition
    !     PGP(:,:,:)  - gridpoint fields (output)

    !                  The ordering of the output fields is as follows (all
    !                  parts are optional depending on the input switches):

    !       vorticity     : IF_UV_G fields
    !       divergence    : IF_UV_G fields
    !       u             : IF_UV_G fields
    !       v             : IF_UV_G fields
    !       scalar fields : IF_SCALARS_G fields
    !       N-S derivative of scalar fields : IF_SCALARS_G fields
    !       E-W derivative of u : IF_UV_G fields
    !       E-W derivative of v : IF_UV_G fields
    !       E-W derivative of scalar fields : IF_SCALARS_G fields

    !     Method.
    !     -------

    !     Externals.  SHUFFLE     - reshuffle fields for load balancing
    !     ----------  FIELD_SPLIT - split fields in NPROMATR packets
    !                 LTINV_CTL   - control of Legendre transform
    !                 FTINV_CTL   - control of Fourier transform

    !     Author.
    !     -------
    !        Mats Hamrud *ECMWF*

    !     Modifications.
    !     --------------
    !        Original : 01-01-03

    !     ------------------------------------------------------------------


    USE PARKIND_ECTRANS,        ONLY: JPIM, JPRB, JPRBT, JPRD
    USE TPM_GEN,                ONLY: NPROMATR
    USE TPM_DISTR,              ONLY: MYSETV
    USE TPM_TRANS,              ONLY: LDIVGP, LSCDERS, LUVDER, LVORGP, GROWING_ALLOCATION, NPROMA, NGPBLKS
    USE ABORT_TRANS_MOD,        ONLY: ABORT_TRANS
    USE BUFFERED_ALLOCATOR_MOD, ONLY: BUFFERED_ALLOCATOR, MAKE_BUFFERED_ALLOCATOR, &
      &                               INSTANTIATE_ALLOCATOR
    USE TRMTOL_MOD,             ONLY: PREPARE_TRMTOL, TRMTOL_HANDLE, TRMTOL
    USE LTINV_VIEW_MOD,         ONLY: PREPARE_LTINV, LTINV_HANDLE, LTINV_VIEW
    USE TRMTOL_PACK_UNPACK,     ONLY: TRMTOL_PACK_HANDLE, TRMTOL_UNPACK_HANDLE, &
      &                               PREPARE_TRMTOL_PACK, PREPARE_TRMTOL_UNPACK, TRMTOL_PACK, &
      &                               TRMTOL_UNPACK
    USE FSC_MOD,                ONLY: FSC
    USE FTINV_MOD,              ONLY: FTINV_HANDLE, PREPARE_FTINV, FTINV
    USE TRLTOG_VIEW_MOD,        ONLY: TRLTOG_HANDLE, PREPARE_TRLTOG, TRLTOG_VIEW
    USE ECTRANS_FIELD_VIEW_INTERNAL_UTIL_MOD, ONLY : SPEC_VIEW, GRID_VIEW

    IMPLICIT NONE

#include "fspgl_intf.h"

    ! Declaration of arguments
    INTEGER(KIND=JPIM) :: KPROMA, KGPBLKS
    TYPE(SPEC_VIEW) :: YDSPVVOR(:), YDSPVDIV(:)
    TYPE(SPEC_VIEW) :: YDSPVSCALAR(:)

    TYPE(GRID_VIEW) :: YDGVU(:),YDGVV(:)
    TYPE(GRID_VIEW) :: YDGVVOR(:),YDGVDIV(:)
    TYPE(GRID_VIEW) :: YDGVSCALAR(:)

    TYPE(GRID_VIEW) :: YDGVU_EW(:),YDGVV_EW(:)
    TYPE(GRID_VIEW) :: YDGVSCALAR_NS(:), YDGVSCALAR_EW(:)

    PROCEDURE (FSPGL_INTF), POINTER, OPTIONAL, INTENT(IN)  :: FSPGL_PROC

    ! Local variables

    REAL(KIND=JPRB), POINTER :: FOUBUF(:), FOUBUF_IN(:)
    REAL(KIND=JPRBT), POINTER :: PREEL_REAL(:), PREEL_COMPLEX(:)
    REAL(KIND=JPRBT), POINTER :: ZOUTS(:), ZOUTA(:)
    REAL(KIND=JPRD), POINTER :: ZOUTS0(:), ZOUTA0(:)
    INTEGER(KIND=JPIM) :: KUV_OFFSET, KSCALARS_OFFSET, KSCALARS_NSDER_OFFSET, &
        & KUV_EWDER_OFFSET, KSCALARS_EWDER_OFFSET
    INTEGER(KIND=JPIM) :: IF_LEG, IF_FOURIER
    INTEGER(KIND=JPIM) :: IFIRST
    TYPE(GRID_VIEW),ALLOCATABLE :: YLGP(:)
    INTEGER(KIND=JPIM) :: IF_UV,IF_UV_G,IF_SCALARS,IF_SCALARS_G,IF_FS,IF_GP,IF_OUT_LT
    INTEGER(KIND=JPIM) :: IF_SCDERS,IF_SCDERS_G
    
    TYPE(BUFFERED_ALLOCATOR) :: ALLOCATOR
    TYPE(LTINV_HANDLE) :: HLTINV
    TYPE(TRMTOL_PACK_HANDLE) :: HTRMTOL_PACK
    TYPE(TRMTOL_HANDLE) :: HTRMTOL
    TYPE(TRMTOL_UNPACK_HANDLE) :: HTRMTOL_UNPACK
    TYPE(FTINV_HANDLE) :: HFTINV
    TYPE(TRLTOG_HANDLE) :: HTRLTOG
    
    INTEGER(KIND=JPIM) :: I,J
    INTEGER(KIND=JPIM), ALLOCATABLE :: IVSET(:)
    !     ------------------------------------------------------------------

    NPROMA = KPROMA
    NGPBLKS = KGPBLKS

    IF (NPROMATR > 0) THEN
      CALL ABORT_TRANS("NPROMATR > 0 not supported for GPU")
    ENDIF
    
    YLGP = [YDGVVOR,YDGVDIV,YDGVU,YDGVV,YDGVSCALAR,YDGVSCALAR_NS,YDGVU_EW,YDGVV_EW,YDGVSCALAR_EW]
    IF_GP = SIZE(YLGP)

    IF_UV_G = SIZE(YDGVU)
    IF_SCALARS_G = SIZE(YDGVSCALAR)
    IF_SCDERS_G = SIZE(YDGVSCALAR_NS)
    ! Compute Vertical domain decomposition

    IF_UV = 0
    DO J = 1, IF_UV_G
      IF (YDGVU(J)%IVSET == MYSETV .OR. YDGVU(J)%IVSET == -1) IF_UV = IF_UV + 1
    ENDDO
    IF_SCALARS = 0
    DO J = 1, IF_SCALARS_G
      IF (YDGVSCALAR(J)%IVSET == MYSETV .OR. YDGVSCALAR(J)%IVSET == -1) IF_SCALARS = IF_SCALARS + 1
    ENDDO
    IF_SCDERS = 0
    DO J = 1, IF_SCDERS_G
      IF (YDGVSCALAR_NS(J)%IVSET == MYSETV .OR. YDGVSCALAR_NS(J)%IVSET == -1) IF_SCDERS = IF_SCDERS + 1
    ENDDO

    ! Initialize potentially unset offsets
    KSCALARS_NSDER_OFFSET = -1
    KUV_EWDER_OFFSET = -1
    KSCALARS_EWDER_OFFSET = -1

    IF_OUT_LT = 2*IF_UV + IF_SCALARS+IF_SCDERS

    LVORGP = SIZE(YDGVVOR) > 0
    LDIVGP = SIZE(YDGVDIV) > 0
    LUVDER = SIZE(YDGVU_EW) > 0 .OR. SIZE(YDGVV_EW) > 0
    LSCDERS = SIZE(YDGVSCALAR_NS) > 0 .OR. SIZE(YDGVSCALAR_EW) > 0

    IF(IF_UV > 0 .AND. LVORGP) THEN
      IF_OUT_LT = IF_OUT_LT+IF_UV
    ENDIF
    IF(IF_UV > 0 .AND. LDIVGP) THEN
      IF_OUT_LT = IF_OUT_LT+IF_UV
    ENDIF
    IF_FS = IF_OUT_LT+IF_SCDERS
    IF(IF_UV > 0 .AND. LUVDER) THEN
      IF_FS = IF_FS+2*IF_UV
    ENDIF

    ! (note in ltinv we will initially start with a slightly different domain decomposition
    ! which always has vorticity and divergence because this is the actual input)
    IFIRST = 0
    IF (LVORGP) IFIRST = IFIRST + IF_UV ! Vorticity
    IF (LDIVGP) IFIRST = IFIRST + IF_UV ! Divergence
    KUV_OFFSET = IFIRST
    IFIRST = IFIRST + IF_UV ! U
    IFIRST = IFIRST + IF_UV ! V
    KSCALARS_OFFSET = IFIRST
    IFIRST = IFIRST + IF_SCALARS ! Scalars
    IF (LSCDERS) THEN
      KSCALARS_NSDER_OFFSET = IFIRST
      IFIRST = IFIRST + IF_SCALARS ! Scalars NS Derivatives
    ENDIF
    ! the rest of fields is being computed  in fourier space, namely in FSC
    IF_LEG = IFIRST
    IF (LUVDER) THEN
      KUV_EWDER_OFFSET = IFIRST
      IFIRST = IFIRST+2*IF_UV ! U and V derivatives
    ENDIF
    IF (LSCDERS) THEN
      KSCALARS_EWDER_OFFSET = IFIRST
      IFIRST = IFIRST + IF_SCALARS ! Scalars EW Derivatives
    ENDIF
    IF_FOURIER = IFIRST
    IF (IF_FOURIER /= IF_FS) CALL ABORT_TRANS('Size mismatch: Wrong computation IF_FS')
    
    ALLOCATOR = MAKE_BUFFERED_ALLOCATOR()
    IF (IF_FS > 0) THEN
      HLTINV = PREPARE_LTINV(ALLOCATOR,IF_UV,IF_SCALARS,LVORGP,LDIVGP,LSCDERS)
      HTRMTOL_PACK = PREPARE_TRMTOL_PACK(ALLOCATOR,IF_LEG)
      HTRMTOL = PREPARE_TRMTOL(ALLOCATOR,IF_LEG)
      HTRMTOL_UNPACK = PREPARE_TRMTOL_UNPACK(ALLOCATOR,IF_FOURIER)
      HFTINV = PREPARE_FTINV(ALLOCATOR,IF_FOURIER)
    ENDIF
    HTRLTOG = PREPARE_TRLTOG(ALLOCATOR,IF_FOURIER,IF_GP)

    CALL INSTANTIATE_ALLOCATOR(ALLOCATOR, GROWING_ALLOCATION)

    IF (IF_FS > 0) THEN
      ! Legendre transformations
      CALL GSTATS(102,0)
      CALL LTINV_VIEW(ALLOCATOR,HLTINV,IF_UV,IF_SCALARS,&
          & YDSPVVOR, YDSPVDIV, YDSPVSCALAR,&
          & ZOUTS,ZOUTA,ZOUTS0,ZOUTA0)
      CALL GSTATS(102,1)

      ! Packing into send buffer, to fourier space and unpack
      CALL GSTATS(152,0)
      CALL TRMTOL_PACK(ALLOCATOR,HTRMTOL_PACK,ZOUTS,ZOUTA,ZOUTS0,ZOUTA0,FOUBUF_IN,IF_LEG)
      CALL TRMTOL(ALLOCATOR,HTRMTOL,FOUBUF_IN,FOUBUF,IF_LEG)
      CALL TRMTOL_UNPACK(ALLOCATOR,HTRMTOL_UNPACK,FOUBUF,PREEL_COMPLEX,IF_LEG,IF_FOURIER)
      CALL GSTATS(152,1)

      CALL GSTATS(107,0)
      ! compute NS derivatives
      CALL FSC(PREEL_COMPLEX, IF_FOURIER, IF_UV, IF_SCALARS, KUV_OFFSET, KSCALARS_OFFSET, &
          & KSCALARS_NSDER_OFFSET, KUV_EWDER_OFFSET, KSCALARS_EWDER_OFFSET)
      !Legendre transformations
      CALL FTINV(ALLOCATOR, HFTINV, PREEL_COMPLEX,PREEL_REAL,IF_FOURIER)
      CALL GSTATS(107,1)
    ENDIF

    ! Transposition into grid-point space
    CALL GSTATS(157,0)
    
    ALLOCATE(IVSET(IF_GP))
    DO I=1,IF_GP
      IVSET(I) = YLGP(I)%IVSET
    ENDDO

    CALL TRLTOG_VIEW(ALLOCATOR,HTRLTOG,PREEL_REAL,IF_FOURIER,IF_GP,IVSET,YLGP)
    DEALLOCATE(IVSET)
    CALL GSTATS(157,1)

  END SUBROUTINE INV_TRANS_VIEW_CTL
END MODULE INV_TRANS_VIEW_CTL_MOD
