SUBROUTINE DIR_TRANSAD(PSPVOR,PSPDIV,PSPSCALAR,PSPSC3A,PSPSC3B,PSPSC2,&
& KPROMA,KVSETUV,KVSETSC,KRESOL,KVSETSC3A,KVSETSC3B,KVSETSC2,&
& PGP,PGPUV,PGP3A,PGP3B,PGP2)


!**** *DIR_TRANSAD* - Direct spectral transform - adjoint.

!     Purpose.
!     --------
!        Interface routine for the direct spectral transform - adjoint

!**   Interface.
!     ----------
!     CALL DIR_TRANSAD(...)

!     Explicit arguments : All arguments except from PGP are optional.
!     -------------------- 
!     PSPVOR(:,:) - spectral vorticity (output)
!     PSPDIV(:,:) - spectral divergence (output)
!     PSPSCALAR(:,:) - spectral scalarvalued fields (output)
!     PSPSC3A(:,:,:) - alternative to use of PSPSCALAR, see PGP3A below (input)
!     PSPSC3B(:,:,:) - alternative to use of PSPSCALAR, see PGP3B below (input)
!     PSPSC2(:,:)  - alternative to use of PSPSCALAR, see PGP2 below (input)
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
!     PGP(:,:,:) - gridpoint fields (input)
!                  PGP need to  dimensioned (NPROMA,IF_GP,NGPBLKS) where
!                  NPROMA is the blocking factor, IF_GP the total number
!                  of output fields and NGPBLKS the number of NPROMA blocks.
!                  The ordering of the output fields is as follows (all 
!                  parts are optional depending on the input switches):
!
!     u             : IF_UV_G fields (if psvor present)
!     v             : IF_UV_G fields (if psvor present)
!     scalar fields : IF_SCALARS_G fields (if pspscalar present)
!   
!     Here IF_UV_G is the GLOBAL number of u/v fields as given by the length
!     of KVSETUV (or by PSPVOR if no split in spectral 'b-set' direction
!     IF_SCALARS_G is the GLOBAL number of scalar fields as giben by the 
!     length of KVESETSC (or by number of fields in PSPSCALAR if no spectral
!     'b-set' split
! 
!     As an alternative to using PGP you can also use a combination of the
!     following arrays. The reason for introducing these alternative ways
!     of calling DIR_TRANS is to avoid uneccessary copies where your data
!     structures don't fit in to the 'PSPVOR,PSPDIV, PSPSCALAR, PGP' layout.
!     The use of any of these precludes the use of PGP and vice versa.

!     PGPUV(:,:,:,:) - the 'u-v' related grid-point variables in the order
!                      described for PGP. The second dimension of PGPUV should
!                      be the same as the "global" first dimension of 
!                      PSPVOR,PSPDIV (in the IFS this is the number of levels)
!                      PGPUV need to be dimensioned(NPROMA,ILEVS,IFLDS,NGPBLKS)
!                      IFLDS is the number of 'variables' (u,v)
!     PGP3A(:,:,:,:) - grid-point array directly connected with PSPSC3A
!                      dimensioned(NPROMA,ILEVS,IFLDS,NGPBLKS)
!                      IFLDS is the number of 'variables' (the same as in
!                      PSPSC3A )
!     PGP3B(:,:,:,:) - grid-point array directly connected with PSPSC3B
!                      dimensioned(NPROMA,ILEVS,IFLDS,NGPBLKS)
!                      IFLDS is the number of 'variables' (the same as in
!                      PSPSC3B)
!     PGP2(:,:,:)    - grid-point array directly connected with PSPSC2
!                      dimensioned(NPROMA,IFLDS,NGPBLKS)
!                      IFLDS is the number of 'variables' (the same as in
!                      PSPSC2 )
! 
!     Method.
!     -------

!     Externals.  SET_RESOL   - set resolution
!     ----------  DIR_TRANS_CTLAD - control routine
!                 

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
REAL_B    ,OPTIONAL, INTENT(INOUT) :: PSPSCALAR(:,:)
REAL_B    ,OPTIONAL, INTENT(INOUT) :: PSPSC3A(:,:,:)
REAL_B    ,OPTIONAL, INTENT(INOUT) :: PSPSC3B(:,:,:)
REAL_B    ,OPTIONAL, INTENT(INOUT) :: PSPSC2(:,:)
INTEGER_M ,OPTIONAL, INTENT(IN) :: KPROMA
INTEGER_M ,OPTIONAL, INTENT(IN) :: KVSETUV(:)
INTEGER_M ,OPTIONAL, INTENT(IN) :: KVSETSC(:)
INTEGER_M ,OPTIONAL, INTENT(IN) :: KVSETSC3A(:)
INTEGER_M ,OPTIONAL, INTENT(IN) :: KVSETSC3B(:)
INTEGER_M ,OPTIONAL, INTENT(IN) :: KVSETSC2(:)
INTEGER_M ,OPTIONAL, INTENT(IN) :: KRESOL

REAL_B,OPTIONAL    ,INTENT(OUT) :: PGP(:,:,:)
REAL_B,OPTIONAL    ,INTENT(OUT) :: PGPUV(:,:,:,:)
REAL_B,OPTIONAL    ,INTENT(OUT) :: PGP3A(:,:,:,:)
REAL_B,OPTIONAL    ,INTENT(OUT) :: PGP3B(:,:,:,:)
REAL_B,OPTIONAL    ,INTENT(OUT) :: PGP2(:,:,:)


END SUBROUTINE DIR_TRANSAD


