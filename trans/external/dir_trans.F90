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
!     Method.
!     -------

!     Externals.  SET_RESOL   - set resolution
!     ----------  LTDIR_CTL   - control of Legendre transform
!                 FTDIR_CTL   - control of Fourier transform

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
USE DIR_TRANS_CTL_MOD
USE ABORT_TRANS_MOD

!endif INTERFACE

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

!ifndef INTERFACE

! Local variables
INTEGER_M :: IUBOUND(3),J
INTEGER_M :: IF_UV,IF_UV_G,IF_SCALARS,IF_SCALARS_G,IF_FS,IF_GP
LOGICAL :: LLPUV,LLPSC

!     ------------------------------------------------------------------

! Set current resolution
CALL SET_RESOL(KRESOL)

! Set defaults

IF_UV = 0
IF_UV_G = 0
IF_SCALARS = 0
IF_SCALARS_G = 0
NPROMA = D%NGPTOT
LLPUV = .FALSE.
LLPSC = .FALSE.
! Decide requirements

IF(PRESENT(KVSETUV)) THEN
  LLPUV = .TRUE.
  IF_UV_G = UBOUND(KVSETUV,1)
  DO J=1,IF_UV_G
    IF(KVSETUV(J) > NPRTRV .OR. KVSETUV(J) < 1) THEN
      WRITE(NERR,*) 'DIR_TRANS:KVSETUV(J) > NPRTRV ',J,KVSETUV(J),NPRTRV
      CALL ABORT_TRANS('DIR_TRANS:KVSETUV TOO LONG OR CONTAINS VALUES OUTSIDE RANGE')
    ENDIF
    IF(KVSETUV(J) == MYSETV) THEN
      IF_UV = IF_UV+1
    ENDIF
  ENDDO
ELSEIF(PRESENT(PSPVOR)) THEN
  IF_UV = UBOUND(PSPVOR,1)
  IF_UV_G = IF_UV
ENDIF

IF(PRESENT(KVSETSC)) THEN
  LLPSC = .TRUE.
  IF_SCALARS_G = UBOUND(KVSETSC,1)
  DO J=1,IF_SCALARS_G
    IF(KVSETSC(J) > NPRTRV .OR. KVSETSC(J) < 1) THEN
      WRITE(NERR,*) 'DIR_TRANS:KVSETSC(J) > NPRTRV ',J,KVSETSC(J),NPRTRV
      CALL ABORT_TRANS('DIR_TRANS:KVSETSC TOO LONG OR CONTAINS VALUES OUTSIDE RANGE')
    ENDIF
    IF(KVSETSC(J) == MYSETV) THEN
      IF_SCALARS = IF_SCALARS+1
    ENDIF
  ENDDO
ELSEIF(PRESENT(PSPSCALAR)) THEN
  IF_SCALARS = UBOUND(PSPSCALAR,1)
  IF_SCALARS_G = IF_SCALARS
ENDIF

IF(PRESENT(KPROMA)) THEN
  NPROMA = KPROMA
ENDIF



! Compute derived variables


NGPBLKS = (D%NGPTOT-1)/NPROMA+1

IF_FS = 2*IF_UV + IF_SCALARS

IF_GP = 2*IF_UV_G+IF_SCALARS_G

! Consistency checks

IF (IF_UV > 0) THEN
  IF(.NOT. PRESENT(PSPVOR) ) THEN
    CALL ABORT_TRANS('DIR_TRANS : IF_UV > 0 BUT PSPVOR MISSING')
  ENDIF
  IF(UBOUND(PSPVOR,1) < IF_UV) THEN
    WRITE(NERR,*)'DIR_TRANS : UBOUND(PSPVOR,1) < IF_UV ',UBOUND(PSPVOR,1),IF_UV
    CALL ABORT_TRANS('DIR_TRANS : PSPVOR TOO SHORT')
  ENDIF
  IF(.NOT. PRESENT(PSPDIV) ) THEN
    CALL ABORT_TRANS('DIR_TRANS : PSPVOR PRESENT BUT PSPDIV MISSING')
  ENDIF
  IF(UBOUND(PSPDIV,1) /= IF_UV) THEN
    WRITE(NERR,*)'DIR_TRANS : UBOUND(PSPDIV,1) < IF_UV ',UBOUND(PSPDIV,1),IF_UV
    CALL ABORT_TRANS('DIR_TRANS : INCONSISTENT FIRST DIM. OF PSPVOR AND PSPDIV')
  ENDIF
ENDIF

IF (IF_SCALARS > 0) THEN
  IF(.NOT. PRESENT(PSPSCALAR) ) THEN
    CALL ABORT_TRANS('DIR_TRANS : IF_SCALARS > 0 BUT PSPSCALAR MISSING')
  ENDIF
  IF(UBOUND(PSPSCALAR,1) < IF_SCALARS) THEN
    WRITE(NERR,*)'DIR_TRANS : UBOUND(PSPSCALAR,1) < IF_SCALARS) ',&
     & UBOUND(PSPSCALAR,1),IF_SCALARS
    CALL ABORT_TRANS('DIR_TRANS : PSPSCALAR TOO SHORT')
  ENDIF
ENDIF

IF(NPRTRV >1) THEN
  IF(IF_UV > 0 .AND. .NOT. PRESENT(KVSETUV)) THEN
    WRITE(NERR,*)'NPRTRV >1 AND IF_UV > 0 AND NOT PRESENT(KVSETUV)',&
                 &NPRTRV,IF_UV
    CALL ABORT_TRANS('DIR_TRANS: SPECIFY VERTICAL SPECTRAL DISTRIBUTION!')
  ENDIF
  IF(IF_SCALARS > 0 .AND. .NOT. PRESENT(KVSETSC)) THEN
    WRITE(NERR,*)'NPRTRV >1 AND IF_SCALARS > 0 AND NOT PRESENT(KVSETSC)',&
                 &NPRTRV,IF_SCALARS
    CALL ABORT_TRANS('DIR_TRANS: SPECIFY VERTICAL SPECTRAL DISTRIBUTION!')
  ENDIF
ENDIF
IUBOUND=UBOUND(PGP)
IF(IUBOUND(1) < NPROMA) THEN
  WRITE(NOUT,*)'DIR_TRANS:FIRST DIM. OF PGP TOO SMALL ',IUBOUND(1),NPROMA
  CALL ABORT_TRANS('DIR_TRANS:FIRST DIMENSION OF PGP TOO SMALL ')
ENDIF
IF(IUBOUND(2) < IF_GP) THEN
  WRITE(NOUT,*)'DIR_TRANS:SEC. DIM. OF PGP TOO SMALL ',IUBOUND(2),IF_GP
  CALL ABORT_TRANS('DIR_TRANS:SECOND DIMENSION OF PGP TOO SMALL ')
ENDIF
IF(IUBOUND(3) < NGPBLKS) THEN
  WRITE(NOUT,*)'DIR_TRANS:THIRD DIM. OF PGP TOO SMALL ',IUBOUND(3),NGPBLKS
  CALL ABORT_TRANS('DIR_TRANS:THIRD DIMENSION OF PGP TOO SMALL ')
ENDIF

!     ------------------------------------------------------------------

CALL DIR_TRANS_CTL(IF_UV_G,IF_SCALARS_G,IF_GP,IF_FS,IF_UV,IF_SCALARS,&
 & PSPVOR,PSPDIV,PSPSCALAR,KVSETUV,KVSETSC,PGP)

!     ------------------------------------------------------------------
!endif INTERFACE

END SUBROUTINE DIR_TRANS

