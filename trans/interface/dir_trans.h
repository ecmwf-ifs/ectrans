SUBROUTINE DIR_TRANS(PSPVOR,PSPDIV,PSPSCALAR,&
& KPROMA,KVSETUV,KVSETSC,KRESOL,&
& PGP)


!**** *DIR_TRANS* - Direct spectral transform (from grid-point to spectral).

!     Purpose.
!     --------
!        Interface routine for the direct spectral transform

!**   Interface.
!     ----------
!     CALL DIR_TRANS(...)

!     Explicit arguments : All arguments except from PGP are optional.
!     -------------------- 
!     PSPVOR(:,:) - spectral vorticity (output)
!     PSPDIV(:,:) - spectral divergence (output)
!     PSPSCALAR(:,:) - spectral scalarvalued fields (output)
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
!     PGP(:,:,:) - gridpoint fields (input)
!                  PGP need to  dimensioned (NPROMA,NF_GP,NGPBLKS) where
!                  NPROMA is the blocking factor, NF_GP the total number
!                  of output fields and NGPBLKS the number of NPROMA blocks.
!                  The ordering of the output fields is as follows (all 
!                  parts are optional depending on the input switches):
!
!       u             : NF_UV_G fields (if psvor present)
!       v             : NF_UV_G fields (if psvor present)
!       scalar fields : NF_SCALARS_G fields (if pspscalar present)
!   
!       Here NF_UV_G is the GLOBAL number of u/v fields as given by the length
!       of KVSETUV (or by PSPVOR if no split in spectral 'b-set' direction
!       NF_SCALARS_G is the GLOBAL number of scalar fields as giben by the 
!       length of KVESETSC (or by number of fields in PSPSCALAR if no spectral
!       'b-set' split
! 
!     ------------------------------------------------------------------

#include "tsmbkind.h"


IMPLICIT NONE

! Declaration of arguments

REAL_B    ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
REAL_B    ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
REAL_B    ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)
INTEGER_M ,OPTIONAL, INTENT(IN) :: KPROMA
INTEGER_M ,OPTIONAL, INTENT(IN) :: KVSETUV(:)
INTEGER_M ,OPTIONAL, INTENT(IN) :: KVSETSC(:)
INTEGER_M ,OPTIONAL, INTENT(IN) :: KRESOL

REAL_B    ,INTENT(IN) :: PGP(:,:,:)


END SUBROUTINE DIR_TRANS

