module suwavedi_mod
contains
SUBROUTINE SUWAVEDI(KSMAX,KTMAX,&
                    &KPROCA,KMYSETA,KASM0,KSPOLEGL,KPROCM,&
                    &KUMPP,KSPEC,KSPEC2,KSPEC2MX,KPOSSP,KMYMS)

!**** *SUWAVEDI * - Routine to initialize spectral wave distribution

!     Purpose.
!     --------
!           Initialize arrays controlling soectral wave distribution

!**   Interface.
!     ----------
!        *CALL* *SUWAVEDI *

!        Explicit arguments : 
!        --------------------
!           KSMAX    - Spectral truncation limit (input)
!           KTMAX    - Overtruncation for KSMAX (input)
!           KPROCA   - Number of processors in A-direction (input)
!           KMYSETA  - A-set for present processor (input)
!           KASM0    - Offsets for spectral waves (output)
!           KSPOLEGL - Local version of NSPOLEG (output)
!           KPROCM   - Where a certain spectral wave belongs  (output)
!           KUMPP    - Number of spectral waves on this PE (output)
!           KSPEC    - Local version on NSPEC (output)
!           KSPEC2   - Local version on NSPEC2 (output)
!           KSPEC2MX - Maximum KSPEC2 across PEs (output)
!           KPOSSP   - Global spectral fields partitioning (output)
!           KMYMS    - This PEs spectral zonal wavenumbers (output)

!        Implicit arguments : NONE
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.   NONE.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        MPP Group *ECMWF*

!     Modifications.
!     --------------
!        Original : 96-01-10
!        L.Isaksen: 96-02-02 - Calculation of KSPEC2MX added
!        K.YESSAD : 97-02-18 - Add KTMAX, bug correction for KSPOLEGL.
!     ------------------------------------------------------------------

#include "tsmbkind.h"

IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KMYSETA
INTEGER_M :: KPROCA
INTEGER_M :: KSMAX
INTEGER_M :: KSPEC
INTEGER_M :: KSPEC2
INTEGER_M :: KSPEC2MX
INTEGER_M :: KSPOLEGL
INTEGER_M :: KTMAX



INTEGER_M :: KASM0(0:KSMAX), KPROCM(0:KSMAX), KUMPP(KPROCA)
INTEGER_M :: KMYMS(KSMAX+1), KPOSSP(KPROCA+1)

INTEGER_M :: ISPEC(KPROCA)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IK, IL, IND, IPOS, ISPEC2P, JA, JM


!      -----------------------------------------------------------------

!*       1.    Initialize partitioning of wave numbers to PEs
!              ----------------------------------------------

ISPEC(:) = 0

KUMPP(:) = 0
KASM0(:)=-99
KSPOLEGL = 0

IL  = 1
IND = 1
IK  = 0
IPOS = 1
DO JM=0,KSMAX
  IK = IK + IND
  IF (IK > KPROCA) THEN
    IK = KPROCA
    IND = -1
  ELSEIF (IK < 1) THEN
    IK = 1
    IND = 1
  ENDIF
  KPROCM(JM) = IK
  ISPEC(IK) = ISPEC(IK)+KSMAX-JM+1
  KUMPP(IK) = KUMPP(IK)+1
  IF (IK == KMYSETA) THEN
    KSPOLEGL = KSPOLEGL +KTMAX+1-JM+1
    KMYMS(IL) = JM
    KASM0(JM) = IPOS
    IPOS = IPOS+(KSMAX-JM+1)*2
    IL = IL+1
  ENDIF
ENDDO

KSPEC  = ISPEC(KMYSETA)
KSPEC2 = 2*KSPEC

KPOSSP(1) = 1
ISPEC2P = 2*ISPEC(1)
KSPEC2MX = ISPEC2P
DO JA=2,KPROCA
  KPOSSP(JA) = KPOSSP(JA-1)+ISPEC2P
  ISPEC2P = 2*ISPEC(JA)
  KSPEC2MX=MAX(KSPEC2MX,ISPEC2P)
ENDDO
KPOSSP(KPROCA+1) = KPOSSP(KPROCA)+ISPEC2P

RETURN
END SUBROUTINE SUWAVEDI
end module suwavedi_mod


