MODULE LDSPC2_MOD
CONTAINS
SUBROUTINE LDSPC2(KM,KF_UV,PEPSNM,POA1,POA2)


USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_DIM

USE UVTVD_MOD


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

IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER(KIND=JPIM), INTENT(IN)  :: KM
INTEGER(KIND=JPIM), INTENT(IN)  :: KF_UV


REAL(KIND=JPRB), INTENT(IN)     :: PEPSNM(0:R%NTMAX+2)
REAL(KIND=JPRB), INTENT(INOUT)  :: POA1(:,:)
REAL(KIND=JPRB), INTENT(OUT)    :: POA2(:,:)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: IUS, IUE, IVS, IVE, IVORS, IVORE, IDIVS, IDIVE


!     ------------------------------------------------------------------



!*       2.   COMPUTE VORTICITY AND DIVERGENCE FROM U AND V.
!             ----------------------------------------------
IUS = 1
IUE = 2*KF_UV
IVS = 2*KF_UV+1
IVE = 4*KF_UV
IVORS = 1
IVORE = 2*KF_UV
IDIVS = 2*KF_UV+1
IDIVE = 4*KF_UV

CALL UVTVD(KM,KF_UV,PEPSNM,POA1(:,IUS:IUE),POA1(:,IVS:IVE),&
          & POA2(:,IVORS:IVORE),POA2(:,IDIVS:IDIVE))

!     ------------------------------------------------------------------

END SUBROUTINE LDSPC2
END MODULE LDSPC2_MOD
