SUBROUTINE INV_TRANS(PSPVOR,PSPDIV,PSPSCALAR,FSPGL_PROC,&
& LDSCDERS,LDVORGP,LDDIVGP,LDUVDER,KPROMA,KVSETUV,KVSETSC,KRESOL,&
& PGP)

!**** *INV_TRANS* - Inverse spectral transform.

!     Purpose.
!     --------
!        Interface routine for the inverse spectral transform

!**   Interface.
!     ----------
!     CALL INV_TRANS(...)

!     Explicit arguments : All arguments except from PGP are optional.
!     -------------------- 
!     PSPVOR(:,:) - spectral vorticity (input)
!     PSPDIV(:,:) - spectral divergence (input)
!     PSPSCALAR(:,:) - spectral scalarvalued fields (input)
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
! 
!     Method.
!     -------

!     Externals.  SET_RESOL   - set resolution
!     ----------  LTINV_CTL   - control of Legendre transform
!                 FTINV_CTL   - control of Fourier transform

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03

!     ------------------------------------------------------------------

#include "tsmbkind.h"


IMPLICIT NONE

! Declaration of arguments

REAL_B    ,OPTIONAL, INTENT(INOUT) :: PSPVOR(:,:)
REAL_B    ,OPTIONAL, INTENT(INOUT) :: PSPDIV(:,:)
REAL_B    ,OPTIONAL, INTENT(IN) :: PSPSCALAR(:,:)
LOGICAL   ,OPTIONAL, INTENT(IN) :: LDSCDERS
LOGICAL   ,OPTIONAL, INTENT(IN) :: LDVORGP
LOGICAL   ,OPTIONAL, INTENT(IN) :: LDDIVGP
LOGICAL   ,OPTIONAL, INTENT(IN) :: LDUVDER
INTEGER_M ,OPTIONAL, INTENT(IN) :: KPROMA
INTEGER_M ,OPTIONAL, INTENT(IN) :: KVSETUV(:)
INTEGER_M ,OPTIONAL, INTENT(IN) :: KVSETSC(:)
INTEGER_M ,OPTIONAL, INTENT(IN) :: KRESOL
EXTERNAL  FSPGL_PROC
OPTIONAL  FSPGL_PROC
REAL_B    ,INTENT(OUT) :: PGP(:,:,:)


END SUBROUTINE INV_TRANS

