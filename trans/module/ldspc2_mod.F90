MODULE LDSPC2_MOD
CONTAINS
SUBROUTINE LDSPC2(KM,PEPSNM,POA1,POA2)


#include "tsmbkind.h"

USE TPM_DIM
USE TPM_TRANS

USE UVTVD_MOD

#ifdef DOC

!**** *LDSPC2* - Spectral computations in direct Legendre transform

!     Purpose.
!     --------
!        To compute vorticity divergence from u and v after the direct
!       Legendre transform.

!**   Interface.
!     ----------
!        CALL LDSPC2(...)

!        Explicit arguments :  
!        --------------------  KM - zonal wave-number
!                              PEPSNM - REPSNM for wavenumber
!                              POA1 - spectral fields for zonal
!                                      wavenumber KM (basic var.)
!                              POA2 - spectral fields for zonal
!                                      wavenumber KM (vor. div.)

!        Implicit arguments :  None
!        --------------------

!     Method.
!     -------

!     Externals.  UVTVD - u and v to vorticity divergence
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 91-07-01
!        Modified : 94-04-06 R. El Khatib - Full-POS configuration 'P'
!        M.Hamrud : 94-11-01 New conf 'G' - vor,div->vor,div
!                            instead of u,v->vor,div
!        MPP Group: 95-10-01 Support for Message Passing version
!     ------------------------------------------------------------------
#endif

IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KM


REAL_B, INTENT(IN)     :: PEPSNM(0:R%NTMAX+2)
REAL_B, INTENT(INOUT)  :: POA1(:,:)
REAL_B, INTENT(OUT)    :: POA2(:,:)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IUS, IUE, IVS, IVE, IVORS, IVORE, IDIVS, IDIVE


!     ------------------------------------------------------------------



!*       2.   COMPUTE VORTICITY AND DIVERGENCE FROM U AND V.
!             ----------------------------------------------
IUS = 1
IUE = 2*NF_UV
IVS = 2*NF_UV+1
IVE = 4*NF_UV
IVORS = 1
IVORE = 2*NF_UV
IDIVS = 2*NF_UV+1
IDIVE = 4*NF_UV

CALL UVTVD(KM,NF_UV,PEPSNM,POA1(:,IUS:IUE),POA1(:,IVS:IVE),&
          & POA2(:,IVORS:IVORE),POA2(:,IDIVS:IDIVE))

!     ------------------------------------------------------------------

END SUBROUTINE LDSPC2
END MODULE LDSPC2_MOD
