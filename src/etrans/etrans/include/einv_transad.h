INTERFACE
SUBROUTINE EINV_TRANSAD(PSPVOR,PSPDIV,PSPSCALAR,PSPSC3A,PSPSC3B,PSPSC2,&
 & FSPGL_PROC,&
 & LDSCDERS,LDVORGP,LDDIVGP,LDUVDER,KPROMA,KVSETUV,KVSETSC,KRESOL,&
 & KVSETSC3A,KVSETSC3B,KVSETSC2,&
 & PGP,PGPUV,PGP3A,PGP3B,PGP2,PMEANU,PMEANV)

!**** *EINV_TRANSAD* - Inverse spectral transform - adjoint.

!     Purpose.
!     --------
!        Interface routine for the inverse spectral transform - adjoint

!**   Interface.
!     ----------
!     CALL EINV_TRANSAD(...)

!     Explicit arguments : All arguments except from PGP are optional.
!     --------------------
!     PSPVOR(:,:) - spectral vorticity (input)
!     PSPDIV(:,:) - spectral divergence (input)
!     PSPSCALAR(:,:) - spectral scalarvalued fields (input)
!     PSPSC3A(:,:,:) - alternative to use of PSPSCALAR, see PGP3A below (input)
!     PSPSC3B(:,:,:) - alternative to use of PSPSCALAR, see PGP3B below (input)
!     PSPSC2(:,:)  - alternative to use of PSPSCALAR, see PGP2 below (input)
!     FSPGL_PROC  - external procedure to be executed in fourier space
!                   before transposition
!     LDSCDERS    - indicating if derivatives of scalar variables are req.
!     LDVORGP     - indicating if grid-point vorticity is req.
!     LDDIVGP     - indicating if grid-point divergence is req.
!     LDUVDER     - indicating if E-W derivatives of u and v are req.
!     KPROMA      - required blocking factor for gridpoint output
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
!     KVSETSC3A(:) - as KVESETSC for PSPSC3A (distribution on first dimension)
!     KVSETSC3B(:) - as KVESETSC for PSPSC3B (distribution on first dimension)
!     KVSETSC2(:) - as KVESETSC for PSPSC2 (distribution on first dimension)
!     KRESOL   - resolution tag  which is required ,default is the
!                first defined resulution (input)
!     PGP(:,:,:) - gridpoint fields (output)
!                  PGP need to  dimensioned (NPROMA,IF_GP,NGPBLKS) where
!                  NPROMA is the blocking factor, IF_GP the total number
!                  of output fields and NGPBLKS the number of NPROMA blocks.
!                  The ordering of the output fields is as follows (all
!                  parts are optional depending on the input switches):
!
!       vorticity     : IF_UV_G fields (if psvor present and LDVORGP)
!       divergence    : IF_UV_G fields (if psvor present and LDDIVGP)
!       u             : IF_UV_G fields (if psvor present)
!       v             : IF_UV_G fields (if psvor present)
!       scalar fields : IF_SCALARS_G fields (if pspscalar present)
!       N-S derivative of scalar fields : IF_SCALARS_G fields (if pspscalar
!                                         present and LDSCDERS)
!       E-W derivative of u : IF_UV_G fields (if psvor present and and LDUVDER)
!       E-W derivative of v : IF_UV_G fields (if psvor present and and LDUVDER)
!       E-W derivative of scalar fields : IF_SCALARS_G fields (if pspscalar
!                                         present and LDSCDERS)
!
!       Here IF_UV_G is the GLOBAL number of u/v fields as given by the length
!       of KVSETUV (or by PSPVOR if no split in spectral 'b-set' direction
!       IF_SCALARS_G is the GLOBAL number of scalar fields as giben by the
!       length of KVESETSC (or by number of fields in PSPSCALAR if no spectral
!       'b-set' split

!     As an alternative to using PGP you can also use a combination of the
!     following arrays. The reason for introducing these alternative ways
!     of calling INV_TRANS is to avoid uneccessary copies where your data
!     structures don't fit in to the 'PSPVOR,PSPDIV, PSPSCALAR, PGP' layout.
!     The use of any of these precludes the use of PGP and vice versa.
!
!     PGPUV(:,:,:,:) - the 'u-v' related grid-point variables in the order
!                      described for PGP. The second dimension of PGPUV should
!                      be the same as the "global" first dimension of
!                      PSPVOR,PSPDIV (in the IFS this is the number of levels)
!                      PGPUV need to be dimensioned(NPROMA,ILEVS,IFLDS,NGPBLKS)
!                      IFLDS is the number of 'variables' (u,v,vor,div ...)
!     PGP3A(:,:,:,:) - grid-point array directly connected with PSPSC3A
!                      dimensioned(NPROMA,ILEVS,IFLDS,NGPBLKS)
!                      IFLDS is the number of 'variables' (the same as in
!                      PSPSC3A if no derivatives, 3 times that with der.)
!     PGP3B(:,:,:,:) - grid-point array directly connected with PSPSC3B
!                      dimensioned(NPROMA,ILEVS,IFLDS,NGPBLKS)
!                      IFLDS is the number of 'variables' (the same as in
!                      PSPSC3B if no derivatives, 3 times that with der.)
!     PGP2(:,:,:)    - grid-point array directly connected with PSPSC2
!                      dimensioned(NPROMA,IFLDS,NGPBLKS)
!                      IFLDS is the number of 'variables' (the same as in
!                      PSPSC2 if no derivatives, 3 times that with der.)

!     Method.
!     -------

!     Externals.  ESET_RESOL   - set resolution
!     ----------  ELTDIR_CTLAD   - control of Legendre transform
!                 EFTDIR_CTLAD   - control of Fourier transform

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03

!     ------------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

! Declaration of arguments

REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSC3A(:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSC3B(:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSC2(:,:)
LOGICAL   ,OPTIONAL, INTENT(IN)  :: LDSCDERS
LOGICAL   ,OPTIONAL, INTENT(IN)  :: LDVORGP
LOGICAL   ,OPTIONAL, INTENT(IN)  :: LDDIVGP
LOGICAL   ,OPTIONAL, INTENT(IN)  :: LDUVDER
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KPROMA
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETUV(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETSC(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETSC3A(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETSC3B(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETSC2(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KRESOL
EXTERNAL  FSPGL_PROC
OPTIONAL  FSPGL_PROC
REAL(KIND=JPRB),OPTIONAL    ,INTENT(IN)  :: PGP(:,:,:)
REAL(KIND=JPRB),OPTIONAL    ,INTENT(IN)  :: PGPUV(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    ,INTENT(IN)  :: PGP3A(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    ,INTENT(IN)  :: PGP3B(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    ,INTENT(IN)  :: PGP2(:,:,:)
REAL(KIND=JPRB),OPTIONAL    ,INTENT(OUT) :: PMEANU(:)
REAL(KIND=JPRB),OPTIONAL    ,INTENT(OUT) :: PMEANV(:)


END SUBROUTINE EINV_TRANSAD

END INTERFACE
