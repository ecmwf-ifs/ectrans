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
!                  PGP need to  dimensioned (NPROMA,NF_GP,NGPBLKS) where
!                  NPROMA is the blocking factor, NF_GP the total number
!                  of output fields and NGPBLKS the number of NPROMA blocks.
!                  The ordering of the output fields is as follows (all 
!                  parts are optional depending on the input switches):
!
!       vorticity     : NF_UV_G fields (if psvor present and LDVORGP)
!       divergence    : NF_UV_G fields (if psvor present and LDDIVGP)
!       u             : NF_UV_G fields (if psvor present)
!       v             : NF_UV_G fields (if psvor present)
!       scalar fields : NF_SCALARS_G fields (if pspscalar present)
!       N-S derivative of scalar fields : NF_SCALARS_G fields (if pspscalar
!                                         present and LDSCDERS)
!       E-W derivative of u : NF_UV_G fields (if psvor present and and LDUVDER)
!       E-W derivative of v : NF_UV_G fields (if psvor present and and LDUVDER)
!       E-W derivative of scalar fields : NF_SCALARS_G fields (if pspscalar
!                                         present and LDSCDERS)
!   
!       Here NF_UV_G is the GLOBAL number of u/v fields as given by the length
!       of KVSETUV (or by PSPVOR if no split in spectral 'b-set' direction
!       NF_SCALARS_G is the GLOBAL number of scalar fields as giben by the 
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

!ifndef INTERFACE

USE TPM_GEN
USE TPM_DIM
USE TPM_TRANS
USE TPM_DISTR
USE TPM_GEOMETRY
USE TPM_FIELDS
USE TPM_FFT

USE SET_RESOL_MOD
USE LTINV_CTL_MOD
USE FTINV_CTL_MOD

!endif INTERFACE

IMPLICIT NONE

! Declaration of arguments

REAL_B    ,OPTIONAL, INTENT(IN) :: PSPVOR(:,:)
REAL_B    ,OPTIONAL, INTENT(IN) :: PSPDIV(:,:)
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

!ifndef INTERFACE

! Local varaibles
INTEGER_M :: IUBOUND(3),J
LOGICAL   :: LLSCDERS
!     ------------------------------------------------------------------

! Set current resolution
CALL SET_RESOL(KRESOL)

! Set defaults

LUV = .FALSE.
LSCALAR = .FALSE.
LVORGP = .FALSE.
LDIVGP = .FALSE.
LUVDER = .FALSE.
NF_UV = 0
NF_SCALARS = 0
NF_SCDERS = 0
NPROMA = D%NGPTOT
LLSCDERS = .FALSE.
! Decide requirements

IF(PRESENT(KVSETUV)) THEN
  NF_UV_G = UBOUND(KVSETUV,1)
  LUV = .TRUE.
  DO J=1,NF_UV_G
    IF(KVSETUV(J) > NPRTRV .OR. KVSETUV(J) < 1) THEN
      WRITE(NERR,*) 'INV_TRANS:KVSETUV(J) > NPRTRV ',J,KVSETUV(J),NPRTRV
      CALL ABOR1('INV_TRANS:KVSETUV TOO LONG OR CONTAINS VALUES OUTSIDE RANGE')
    ENDIF
    IF(KVSETUV(J) == MYSETV) THEN
      NF_UV = NF_UV+1
    ENDIF
  ENDDO
ELSEIF(PRESENT(PSPVOR)) THEN
  LUV = .TRUE.
  NF_UV = UBOUND(PSPVOR,1)
  NF_UV_G = NF_UV
ENDIF

IF(PRESENT(KVSETSC)) THEN
  LSCALAR = .TRUE.
  NF_SCALARS_G = UBOUND(KVSETSC,1)
  DO J=1,NF_SCALARS_G
    IF(KVSETSC(J) > NPRTRV .OR. KVSETSC(J) < 1) THEN
      WRITE(NERR,*) 'INV_TRANS:KVSETSC(J) > NPRTRV ',J,KVSETSC(J),NPRTRV
      CALL ABOR1('INV_TRANS:KVSETSC TOO LONG OR CONTAINS VALUES OUTSIDE RANGE')
    ENDIF
    IF(KVSETSC(J) == MYSETV) THEN
      NF_SCALARS = NF_SCALARS+1
    ENDIF
  ENDDO
ELSEIF(PRESENT(PSPSCALAR)) THEN
  LSCALAR = .TRUE.
  NF_SCALARS = UBOUND(PSPSCALAR,1)
  NF_SCALARS_G = NF_SCALARS
ENDIF

IF (LSCALAR) THEN
  IF(PRESENT(LDSCDERS)) THEN
    NF_SCDERS = NF_SCALARS
    LLSCDERS = .TRUE.
  ENDIF
ENDIF

IF(PRESENT(KPROMA)) THEN
  NPROMA = KPROMA
ENDIF

IF(PRESENT(LDVORGP)) THEN
  LVORGP = LDVORGP
ENDIF

IF(PRESENT(LDDIVGP)) THEN
  LDIVGP = LDDIVGP
ENDIF

IF(PRESENT(LDUVDER)) THEN
  LUVDER = LDUVDER
ENDIF



! Compute derived variables


IF(LVORGP) LDIVGP = .TRUE.

NGPBLKS = (D%NGPTOT-1)/NPROMA+1

NF_OUT_LT = 2*NF_UV + NF_SCALARS+NF_SCDERS

IF(NF_UV > 0 .AND. LVORGP) THEN
  NF_OUT_LT = NF_OUT_LT+NF_UV
ENDIF
IF(NF_UV > 0 .AND. LDIVGP) THEN
  NF_OUT_LT = NF_OUT_LT+NF_UV
ENDIF
NF_FS = NF_OUT_LT+NF_SCDERS
IF(NF_UV > 0 .AND. LUVDER) THEN
  NF_FS = NF_FS+2*NF_UV
ENDIF

NF_GP = 2*NF_UV_G+NF_SCALARS_G
IF(LLSCDERS) THEN
  NF_GP  = NF_GP+2*NF_SCALARS_G
ENDIF
IF(NF_UV_G > 0 .AND. LVORGP) THEN
  NF_GP = NF_GP+NF_UV_G
ENDIF
IF(NF_UV_G > 0 .AND. LDIVGP) THEN
  NF_GP = NF_GP+NF_UV_G
ENDIF
IF(NF_UV_G > 0 .AND. LUVDER) THEN
  NF_GP = NF_GP+2*NF_UV_G
ENDIF

! Consistency checks

IF (NF_UV > 0) THEN
  IF(.NOT. PRESENT(PSPVOR) ) THEN
    CALL ABOR1("INV_TRANS : NF_UV > 0 BUT PSPVOR MISSING")
  ENDIF
  IF(UBOUND(PSPVOR,1) < NF_UV) THEN
    CALL ABOR1("INV_TRANS : PSPVOR TOO SHORT")
  ENDIF
  IF(.NOT. PRESENT(PSPDIV) ) THEN
    CALL ABOR1("INV_TRANS : NF_UV > 0 BUT PSPDIV MISSING")
  ENDIF
  IF(UBOUND(PSPDIV,1) < NF_UV) THEN
    CALL ABOR1("INV_TRANS : PSPDIV TOO SHORT")
  ENDIF
ENDIF

IF(NF_UV_G == 0) THEN
  LUVDER = .FALSE.
ENDIF
IF(NPRTRV >1) THEN
  IF(NF_UV > 0 .AND. .NOT. PRESENT(KVSETUV)) THEN
    WRITE(NERR,*)'NPRTRV >1 AND NF_UV > 0 AND NOT PRESENT(KVSETUV)',&
                 &NPRTRV,NF_UV
    CALL ABOR1('INV_TRANS: SPECIFY VERTICAL SPECTRAL DISTRIBUTION!')
  ENDIF
  IF(NF_SCALARS > 0 .AND. .NOT. PRESENT(KVSETSC)) THEN
    WRITE(NERR,*)'NPRTRV >1 AND NF_SCALARS > 0 AND NOT PRESENT(KVSETSC)',&
                 &NPRTRV,NF_SCALARS
    CALL ABOR1('INV_TRANS: SPECIFY VERTICAL SPECTRAL DISTRIBUTION!')
  ENDIF
ENDIF
IUBOUND=UBOUND(PGP)
IF(IUBOUND(1) < NPROMA) THEN
  WRITE(NOUT,*)'INV_TRANS:FIRST DIM. OF PGP TOO SMALL ',IUBOUND(1),NPROMA
  CALL ABOR1('INV_TRANS:FIRST DIMENSION OF PGP TOO SMALL ')
ENDIF
IF(IUBOUND(2) < NF_GP) THEN
  WRITE(NOUT,*)'INV_TRANS:SEC. DIM. OF PGP TOO SMALL ',IUBOUND(2),NF_GP
  CALL ABOR1('INV_TRANS:SECOND DIMENSION OF PGP TOO SMALL ')
ENDIF
IF(IUBOUND(3) < NGPBLKS) THEN
  WRITE(NOUT,*)'INV_TRANS:THIRD DIM. OF PGP TOO SMALL ',IUBOUND(3),NGPBLKS
  CALL ABOR1('INV_TRANS:THIRD DIMENSION OF PGP TOO SMALL ')
ENDIF

!     ------------------------------------------------------------------

! Perform transform

IF (NF_OUT_LT > 0) THEN
  CALL LTINV_CTL(PSPVOR,PSPDIV,PSPSCALAR,FSPGL_PROC)
ENDIF

CALL FTINV_CTL(PGP,KVSETUV,KVSETSC)

!     ------------------------------------------------------------------

!endif INTERFACE

END SUBROUTINE INV_TRANS

