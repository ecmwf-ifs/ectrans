MODULE LTDIR_MOD
CONTAINS
SUBROUTINE LTDIR(KM,KMLOC,PSPVOR,PSPDIV,PSPSCALAR)


#include "tsmbkind.h"

USE TPM_DIM
USE TPM_TRANS

USE PRLE2_MOD
USE PREPSNM_MOD
USE PRFI2_MOD
USE LDFOU2_MOD
USE LEDIR_MOD
USE LDSPC2_MOD
USE UPDSP_MOD  
 
#ifdef DOC

!**** *LTDIR* - Control of Direct Legendre transform step

!     Purpose.
!     --------
!        Tranform from Fourier space to spectral space, compute
!        vorticity and divergence.

!**   Interface.
!     ----------
!        *CALL* *LTDIR(...)*

!        Explicit arguments : 
!        --------------------  KM     - zonal wavenumber
!                              KMLOC  - local zonal wavenumber

!        Implicit arguments :  None
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------
!         PRLE2   - prepares the Legendre polonymials for truncation NTMAX.
!         PREPSNM - prepare REPSNM for wavenumber KM
!         PRFI2   - prepares the Fourier work arrays for model variables.
!         LDFOU2  - computations in Fourier space
!         LEDIR   - direct Legendre transform
!         LDSPC2  - computations in spectral space
!         UPDSP   - updating of spectral arrays (fields)

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 87-11-24
!        Modified : 91-07-01 Philippe Courtier/Mats Hamrud - Rewrite
!                            for uv formulation
!        Modified 93-03-19 D. Giard - CDCONF='T' for tendencies
!        Modified 93-11-18 M. Hamrud - use only one Fourier buffer
!        Modified 94-04-06 R. El khatib Full-POS implementation
!        M.Hamrud  : 94-11-01 New conf 'G' - vor,div->vor,div
!                             instead of u,v->vor,div
!        MPP Group : 95-10-01 Support for Distributed Memory version
!        K. YESSAD (AUGUST 1996): 
!               - Legendre transforms for transmission coefficients.
!        Modified : 04/06/99 D.Salmond : change order of AIA and SIA
!     ------------------------------------------------------------------
#endif

IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M, INTENT(IN)  :: KM
INTEGER_M, INTENT(IN)  :: KMLOC

REAL_B  ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
REAL_B  ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
REAL_B  ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)


!     LOCAL INTEGER SCALARS
INTEGER_M :: IFC

!     LOCAL REALS
REAL_B :: ZSIA(NLED2,R%NDGNH),       ZAIA(NLED2,R%NDGNH)
REAL_B :: ZLEPO(R%NLED3,R%NDGNH)
REAL_B :: ZEPSNM(0:R%NTMAX+2)
REAL_B :: ZOA1(R%NLED4,NLED2),         ZOA2(R%NLED4,NLED2)


!     ------------------------------------------------------------------

!*       1.    PREPARE LEGENDRE POLONOMIALS AND EPSNM
!              --------------------------------------


CALL PRLE2(KM,ZLEPO)

CALL PREPSNM(KM,KMLOC,ZEPSNM)

!     ------------------------------------------------------------------

!*       2.    PREPARE WORK ARRAYS.
!              --------------------

CALL PRFI2(KM,KMLOC,ZAIA,ZSIA)

!     ------------------------------------------------------------------

!*       3.    FOURIER SPACE COMPUTATIONS.
!              ---------------------------

CALL LDFOU2(KM,ZAIA,ZSIA)

!     ------------------------------------------------------------------

!*       4.    DIRECT LEGENDRE TRANSFORM.
!              --------------------------
IFC = 2*NF_FS 
CALL LEDIR(KM,IFC,ZAIA,ZSIA,ZOA1,ZLEPO)

!     ------------------------------------------------------------------

!*       5.    COMPUTE VORTICITY AND DIVERGENCE.
!              ---------------------------------

CALL LDSPC2(KM,ZEPSNM,ZOA1,ZOA2)

!     ------------------------------------------------------------------

!*       6.    UPDATE SPECTRAL ARRAYS.
!              -----------------------

CALL UPDSP(KM,ZOA1,ZOA2,PSPVOR,PSPDIV,PSPSCALAR)

!     ------------------------------------------------------------------

END SUBROUTINE LTDIR
END MODULE LTDIR_MOD

