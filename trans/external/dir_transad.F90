SUBROUTINE DIR_TRANSAD(PSPVOR,PSPDIV,PSPSCALAR,&
& KPROMA,KVSETUV,KVSETSC,KRESOL,&
& PGP)


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
!     u             : NF_UV_G fields (if psvor present)
!     v             : NF_UV_G fields (if psvor present)
!     scalar fields : NF_SCALARS_G fields (if pspscalar present)
!   
!     Here NF_UV_G is the GLOBAL number of u/v fields as given by the length
!     of KVSETUV (or by PSPVOR if no split in spectral 'b-set' direction
!     NF_SCALARS_G is the GLOBAL number of scalar fields as giben by the 
!     length of KVESETSC (or by number of fields in PSPSCALAR if no spectral
!     'b-set' split
! 
!     Method.
!     -------

!     Externals.  SET_RESOL   - set resolution
!     ----------  LTDIR_CTLAD - control of Legendre transform
!                 FTDIR_CTLAD - control of Fourier transform

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03

!     ------------------------------------------------------------------

#include "tsmbkind.h"

!ifndef INTERFACE

USE TPM_GEN
USE TPM_TRANS
USE TPM_DISTR

USE SET_RESOL_MOD
USE LTDIR_CTLAD_MOD
USE FTDIR_CTLAD_MOD

!endif INTERFACE

IMPLICIT NONE

! Declaration of arguments

REAL_B    ,OPTIONAL, INTENT(INOUT) :: PSPVOR(:,:)
REAL_B    ,OPTIONAL, INTENT(INOUT) :: PSPDIV(:,:)
REAL_B    ,OPTIONAL, INTENT(INOUT) :: PSPSCALAR(:,:)
INTEGER_M ,OPTIONAL, INTENT(IN) :: KPROMA
INTEGER_M ,OPTIONAL, INTENT(IN) :: KVSETUV(:)
INTEGER_M ,OPTIONAL, INTENT(IN) :: KVSETSC(:)
INTEGER_M ,OPTIONAL, INTENT(IN) :: KRESOL

REAL_B    ,INTENT(OUT) :: PGP(:,:,:)

!ifndef INTERFACE

! Local variables
INTEGER_M :: IUBOUND(3),J

!     ------------------------------------------------------------------

! Set current resolution
CALL SET_RESOL(KRESOL)

! Set defaults

LUV = .FALSE.
LSCALAR = .FALSE.
NF_UV = 0
NF_SCALARS = 0
NPROMA = D%NGPTOT
! Decide requirements

IF(PRESENT(PSPVOR)) THEN
  LUV = .TRUE.
  NF_UV = UBOUND(PSPVOR,1)
ENDIF

IF(PRESENT(PSPSCALAR)) THEN
  LSCALAR = .TRUE.
  NF_SCALARS = UBOUND(PSPSCALAR,1)
ENDIF

IF(PRESENT(KPROMA)) THEN
  NPROMA = KPROMA
ENDIF


IF(PRESENT(KVSETUV)) THEN
  NF_UV_G = UBOUND(KVSETUV,1)
  DO J=1,NF_UV_G
    IF(KVSETUV(J) > NPRTRV) THEN
      WRITE(NERR,*) 'DIR_TRANSAD:KVSETUV(J) > NPRTRV ',J,KVSETUV(J),NPRTRV
      CALL ABOR1('DIR_TRANSAD:KVSETUV  CONTAINS VALUES OUTSIDE RANGE')
    ENDIF
  ENDDO
ELSE
  NF_UV_G = NF_UV
ENDIF

IF(PRESENT(KVSETSC)) THEN
  NF_SCALARS_G = UBOUND(KVSETSC,1)
  DO J=1,NF_SCALARS_G
    IF(KVSETSC(J) > NPRTRV) THEN
      WRITE(NERR,*) 'DIR_TRANSAD:KVSETSC(J) > NPRTRV ',J,KVSETSC(J),NPRTRV
      CALL ABOR1('DIR_TRANSAD:KVSETSC CONTAINS VALUES OUTSIDE RANGE')
    ENDIF
  ENDDO
ELSE
  NF_SCALARS_G = NF_SCALARS
ENDIF

! Compute derived variables


NGPBLKS = (D%NGPTOT-1)/NPROMA+1

NF_FS = 2*NF_UV + NF_SCALARS

NF_GP = 2*NF_UV_G+NF_SCALARS_G

! Consistency checks

IF (LUV) THEN
  IF(.NOT. PRESENT(PSPDIV) ) THEN
    CALL ABOR1("DIR_TRANSAD : PSPVOR PRESENT BUT PSPDIV MISSING")
  ENDIF
  IF(UBOUND(PSPDIV,1) /= NF_UV) THEN
    CALL ABOR1("DIR_TRANSAD : INCONSISTENT FIRST DIM. OF PSPVOR AND PSPDIV")
  ENDIF
ENDIF
IF(NF_UV_G < NF_UV ) THEN
  WRITE(NERR,*)'DIR_TRANSAD:NF_UV_G < NF_UV',NF_UV_G,NF_UV
  CALL ABOR1('DIR_TRANSAD:NF_UV_G < NF_UV')
ENDIF
IF(NF_SCALARS_G < NF_SCALARS ) THEN
  WRITE(NERR,*)'DIR_TRANSAD:NF_SCALARS_G < NF_SCALARS',NF_SCALARS_G,NF_SCALARS
  CALL ABOR1('DIR_TRANSAD:NF_SCALARS_G < NF_SCALARS')
ENDIF

IF(NPRTRV >1) THEN
  IF(NF_UV > 0 .AND. .NOT. PRESENT(KVSETUV)) THEN
    WRITE(NERR,*)'NPRTRV >1 AND NF_UV > 0 AND NOT PRESENT(KVSETUV)',&
                 &NPRTRV,NF_UV
    CALL ABOR1('DIR_TRANSAD: SPECIFY VERTICAL SPECTRAL DISTRIBUTION!')
  ENDIF
  IF(NF_SCALARS > 0 .AND. .NOT. PRESENT(KVSETSC)) THEN
    WRITE(NERR,*)'NPRTRV >1 AND NF_SCALARS > 0 AND NOT PRESENT(KVSETSC)',&
                 &NPRTRV,NF_SCALARS
    CALL ABOR1('DIR_TRANSAD: SPECIFY VERTICAL SPECTRAL DISTRIBUTION!')
  ENDIF
ENDIF
IUBOUND=UBOUND(PGP)
IF(IUBOUND(1) < NPROMA) THEN
  WRITE(NOUT,*)'DIR_TRANSAD:FIRST DIM. OF PGP TOO SMALL ',IUBOUND(1),NPROMA
  CALL ABOR1('DIR_TRANSAD:FIRST DIMENSION OF PGP TOO SMALL ')
ENDIF
IF(IUBOUND(2) < NF_GP) THEN
  WRITE(NOUT,*)'DIR_TRANSAD:SEC. DIM. OF PGP TOO SMALL ',IUBOUND(2),NF_GP
  CALL ABOR1('DIR_TRANSAD:SECOND DIMENSION OF PGP TOO SMALL ')
ENDIF
IF(IUBOUND(3) < NGPBLKS) THEN
  WRITE(NOUT,*)'DIR_TRANSAD:THIRD DIM. OF PGP TOO SMALL ',IUBOUND(3),NGPBLKS
  CALL ABOR1('DIR_TRANSAD:THIRD DIMENSION OF PGP TOO SMALL ')
ENDIF

!     ------------------------------------------------------------------

! Perform transform

CALL LTDIR_CTLAD(PSPVOR,PSPDIV,PSPSCALAR)

CALL FTDIR_CTLAD(PGP,KVSETUV,KVSETSC)


!     ------------------------------------------------------------------
!endif INTERFACE

END SUBROUTINE DIR_TRANSAD


