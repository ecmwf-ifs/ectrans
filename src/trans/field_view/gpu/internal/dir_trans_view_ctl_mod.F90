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

MODULE DIR_TRANS_VIEW_CTL_MOD
CONTAINS
  SUBROUTINE DIR_TRANS_VIEW_CTL(KPROMA,KGPBLKS, &
                              & YDSPVVOR, YDSPVDIV, YDSPVSCALAR, &
                              & YDGVU,YDGVV,YDGVSCALAR)

    !**** *DIR_TRANS_VIEW_CTL* - Control routine for direct spectral transform.

    !     Purpose.
    !     --------
    !        Control routine for the direct spectral transform

    !**   Interface.
    !     ----------
    !     CALL DIR_TRANS_VIEW_CTL(...)

    !     Explicit arguments :
    !     --------------------
    !     KPROMA       - global number of spectral u-v fields
    !     KGPBLKS      - global number of scalar spectral fields    
    !     IF_FS        - total number of fields in fourier space
    !     KF_UV        - local number of spectral u-v fields
    !     KF_SCALARS   - local number of scalar spectral fields
    !     PSPVOR(:,:)  - spectral vorticity
    !     PSPDIV(:,:)  - spectral divergence
    !     PSPSCALAR(:,:) - spectral scalarvalued fields
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
    !     PGP(:,:,:)  - gridpoint fields

    !                  The ordering of the output fields is as follows (all
    !                  parts are optional depending on the input switches):
    !
    !       u             : KF_UV_G fields
    !       v             : KF_UV_G fields
    !       scalar fields : KF_SCALARS_G fields

    !     Method.
    !     -------

    !     Externals.  SHUFFLE     - reshuffle fields for load balancing
    !     ----------  FIELD_SPLIT - split fields in NPROMATR packets
    !                 LTDIR_CTL   - control of Legendre transform
    !                 FTDIR_CTL   - control of Fourier transform

    !     Author.
    !     -------
    !        Mats Hamrud *ECMWF*

    !     Modifications.
    !     --------------
    !        Original : 01-01-03

    !     ------------------------------------------------------------------

    USE PARKIND_ECTRANS,        ONLY: JPRBT, JPRD, JPRB, JPIM
    USE TPM_GEN,                ONLY: NPROMATR
    USE TPM_TRANS,              ONLY: GROWING_ALLOCATION
    USE BUFFERED_ALLOCATOR_MOD, ONLY: BUFFERED_ALLOCATOR, MAKE_BUFFERED_ALLOCATOR, &
      &                               INSTANTIATE_ALLOCATOR
    USE FTDIR_MOD,              ONLY: FTDIR_HANDLE, PREPARE_FTDIR, FTDIR
    USE LTDIR_VIEW_MOD,         ONLY: LTDIR_HANDLE, PREPARE_LTDIR, LTDIR_VIEW
    USE TRGTOL_VIEW_MOD,        ONLY: TRGTOL_HANDLE, PREPARE_TRGTOL, TRGTOL_VIEW
    USE TRLTOM_MOD,             ONLY: TRLTOM_HANDLE, PREPARE_TRLTOM, TRLTOM
    USE TRLTOM_PACK_UNPACK,     ONLY: TRLTOM_PACK_HANDLE, TRLTOM_UNPACK_HANDLE, &
      &                               PREPARE_TRLTOM_PACK, PREPARE_TRLTOM_UNPACK, TRLTOM_PACK, &
      &                               TRLTOM_UNPACK
    USE ABORT_TRANS_MOD,        ONLY: ABORT_TRANS
    USE ECTRANS_FIELD_VIEW_INTERNAL_UTIL_MOD, ONLY : SPEC_VIEW, GRID_VIEW
    USE TPM_TRANS       ,ONLY : NPROMA, NGPBLKS
    USE TPM_DISTR       ,ONLY : MYSETV

    IMPLICIT NONE

    ! Declaration of arguments

    ! Declaration of arguments
    INTEGER(KIND=JPIM) :: KPROMA, KGPBLKS
    TYPE(SPEC_VIEW) :: YDSPVVOR(:), YDSPVDIV(:)
    TYPE(SPEC_VIEW) :: YDSPVSCALAR(:)

    TYPE(GRID_VIEW) :: YDGVU(:),YDGVV(:)
    TYPE(GRID_VIEW) :: YDGVSCALAR(:)
    
    ! Local variables

    !New one    
    INTEGER(KIND=JPIM) :: IF_FS,IF_GP, IF_UV, IF_SCALAR, J
    TYPE(GRID_VIEW),ALLOCATABLE :: YLGP(:)
    INTEGER(KIND=JPIM), ALLOCATABLE :: IVSET(:)

    REAL(KIND=JPRBT), POINTER :: FOUBUF_IN(:), FOUBUF(:)
    REAL(KIND=JPRBT), POINTER :: PREEL_REAL(:), PREEL_COMPLEX(:)

    REAL(KIND=JPRBT), POINTER :: ZINPS(:), ZINPA(:)
    REAL(KIND=JPRD), POINTER :: ZINPS0(:), ZINPA0(:)

    TYPE(BUFFERED_ALLOCATOR) :: ALLOCATOR
    TYPE(TRGTOL_HANDLE) :: HTRGTOL
    TYPE(FTDIR_HANDLE) :: HFTDIR
    TYPE(TRLTOM_PACK_HANDLE) :: HTRLTOM_PACK
    TYPE(TRLTOM_HANDLE) :: HTRLTOM
    TYPE(TRLTOM_UNPACK_HANDLE) :: HTRLTOM_UNPACK
    TYPE(LTDIR_HANDLE) :: HLTDIR

    ! Perform transform
    NPROMA = KPROMA
    NGPBLKS = KGPBLKS

    IF (NPROMATR > 0) THEN
      CALL ABORT_TRANS("NPROMATR > 0 not supported for GPU")
    ENDIF

    YLGP = [YDGVU, YDGVV,YDGVSCALAR]    
    IF_GP = SIZE(YLGP)
    IF_UV = SIZE(YDSPVVOR)
    IF_SCALAR = SIZE(YDSPVSCALAR)

    YLGP = [YDGVU, YDGVV,YDGVSCALAR]
    IF_GP = SIZE(YLGP)
    
      
    ALLOCATE(IVSET(IF_GP))    
    IF_FS = 0
    DO J = 1, IF_GP
      IF (YLGP(J)%IVSET == MYSETV .OR. YLGP(J)%IVSET == -1) IF_FS = IF_FS + 1
      IVSET(J) = YLGP(J)%IVSET
    ENDDO
        
    ! Prepare everything
    ALLOCATOR = MAKE_BUFFERED_ALLOCATOR()
    HTRGTOL = PREPARE_TRGTOL(ALLOCATOR,IF_GP,IF_FS)
    IF (IF_FS > 0) THEN
      HFTDIR = PREPARE_FTDIR(ALLOCATOR,IF_FS)
      HTRLTOM_PACK = PREPARE_TRLTOM_PACK(ALLOCATOR, IF_FS)
      HTRLTOM = PREPARE_TRLTOM(ALLOCATOR, IF_FS)
      HTRLTOM_UNPACK = PREPARE_TRLTOM_UNPACK(ALLOCATOR, IF_FS)
      HLTDIR = PREPARE_LTDIR(ALLOCATOR, IF_FS, IF_UV)
    ENDIF

    CALL INSTANTIATE_ALLOCATOR(ALLOCATOR, GROWING_ALLOCATION)
  
    ! from the PGP arrays to PREEL_REAL
    CALL GSTATS(158,0)
    CALL TRGTOL_VIEW(ALLOCATOR,HTRGTOL,PREEL_REAL,IF_FS,IF_GP,&
     & KVSET = IVSET, YDGP = YLGP)
     
    CALL GSTATS(158,1)

    IF (IF_FS > 0) THEN

      ! fourier transform from PREEL_REAL to PREEL_COMPLEX (in-place!)
      CALL GSTATS(106,0)
      CALL FTDIR(ALLOCATOR,HFTDIR,PREEL_REAL,PREEL_COMPLEX,IF_FS)
      CALL GSTATS(106,1)

      CALL GSTATS(153,0)

      CALL TRLTOM_PACK(ALLOCATOR,HTRLTOM_PACK,PREEL_COMPLEX,FOUBUF_IN,IF_FS)
      CALL TRLTOM(ALLOCATOR,HTRLTOM,FOUBUF_IN,FOUBUF,IF_FS)
      CALL TRLTOM_UNPACK(ALLOCATOR,HTRLTOM_UNPACK,FOUBUF,ZINPS,ZINPA,ZINPS0,ZINPA0,IF_FS,IF_UV)
      CALL GSTATS(153,1)

      CALL GSTATS(103,0)
      CALL LTDIR_VIEW(ALLOCATOR,HLTDIR,ZINPS,ZINPA,ZINPS0,ZINPA0,IF_GP,IF_UV,IF_SCALAR,YDSPVVOR,YDSPVDIV,YDSPVSCALAR)
      CALL GSTATS(103,1)

    ENDIF

    DEALLOCATE(IVSET)

  END SUBROUTINE DIR_TRANS_VIEW_CTL
END MODULE DIR_TRANS_VIEW_CTL_MOD
