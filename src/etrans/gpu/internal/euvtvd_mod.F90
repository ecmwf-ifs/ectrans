MODULE EUVTVD_MOD
CONTAINS
SUBROUTINE EUVTVD(KFIELD,PU,PV,PVOR,PDIV)

!**** *EUVTVD* - Compute vor/div from u and v in spectral space

!     Purpose.
!     --------
!        To compute vorticity and divergence from u and v in spectral
!       space. Input u and v from KM to NTMAX+1, output vorticity and
!       divergence from KM to NTMAX - calculation part.

!**   Interface.
!     ----------
!        CALL EUVTVD(KM,KFIELD,PEPSNM,PU,PV,PVOR,PDIV)

!        Explicit arguments :  KM - zonal wave-number
!        --------------------  KFIELD - number of fields (levels)
!                              KFLDPTR - fields pointers
!                              PEPSNM - REPSNM for wavenumber KM
!                              PU - u wind component for zonal
!                                   wavenumber KM
!                              PV - v wind component for zonal
!                                   wavenumber KM
!                              PVOR - vorticity for zonal
!                                     wavenumber KM
!                              PDIV - divergence for zonal
!                                     wavenumber KM

!     Method.  See ref.
!     -------

!     Externals.  None.
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
!        D. Giard : NTMAX instead of NSMAX
!        01-08-27 : R. El Khatib Fix for NPROMATR /= 0
!        03-03-03 : G. Radnoti: b-level conform mean-wind distribution
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        F. Vana + NEC 28-Apr-2009 MPI-OpenMP fix
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        N. Lopes & R. El Khatib 15-Jun-2012 Scalability enhancement
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DIM         ,ONLY : R
USE TPM_DISTR       ,ONLY : D, NPRCIDS, NPRTRW, MYSETV, MYSETW, MYPROC, NPROC
USE TPM_DISTR       ,ONLY : D_NUMP,D_MYMS
USE TPMALD_GEO      ,ONLY : GALD
USE TPMALD_DISTR    ,ONLY : DALD, DALD_NCPL2M
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
!

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KFIELD
REAL(KIND=JPRB), INTENT(INOUT) :: PU  (:,:,:),PV  (:,:,:)
REAL(KIND=JPRB), INTENT(OUT)   :: PVOR(:,:,:),PDIV(:,:,:)

INTEGER(KIND=JPIM) :: II, IN, IR, J, JN
INTEGER(KIND=JPIM) :: IM, JM, JNMAX

REAL(KIND=JPRB) :: ZKM
REAL(KIND=JPRB) :: ZIN
INTEGER(KIND=JPIM) :: JA,ITAG,ILEN,IFLD,ISND
character(len=64) :: frmt
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EUVTVD_MOD:EUVTVD',0,ZHOOK_HANDLE)

!*       1.    COMPUTE U V FROM VORTICITY AND DIVERGENCE.
!              ------------------------------------------


!$acc parallel loop collapse(3) private(J,JM,JN,IR,II,IM,ZKM) present (PVOR, PDIV, PU, PV,D_MYMS,D_NUMP)
DO J=1,KFIELD
  DO JM=1,D_NUMP
    DO JN=1,R%NDGL+R%NNOEXTZG
      IM = D_MYMS(JM)
      ZKM=REAL(IM,JPRB)*GALD%EXWN
      IR=2*J-1
      II=IR+1
      PDIV(JN,JM,IR)=-ZKM*PU(JN,JM,II)
      PDIV(JN,JM,II)= ZKM*PU(JN,JM,IR)
      PVOR(JN,JM,IR)=-ZKM*PV(JN,JM,II)
      PVOR(JN,JM,II)= ZKM*PV(JN,JM,IR)
    ENDDO
  ENDDO
ENDDO
!$acc end parallel loop

JNMAX = MAXVAL(DALD%NCPL2M)

!$acc parallel loop collapse(3) private(J,JM,JN,IM,ZIN,IN) copyin (JNMAX) present (PVOR, PDIV, PU, PV,DALD_NCPL2M)
DO J=1,2*KFIELD
  DO JM=1,D_NUMP
    DO JN=1,JNMAX,2
      IM = D_MYMS(JM)
!      IF ( JN <= DALD_NCPL2M(IM) ) THEN
        IN=(JN-1)/2
        ZIN=REAL(IN,JPRB)*GALD%EYWN
        PVOR(JN  ,JM,J)=PVOR(JN  ,JM,J)+ZIN*PU(JN+1,JM,J)
        PVOR(JN+1,JM,J)=PVOR(JN+1,JM,J)-ZIN*PU(JN  ,JM,J)
        PDIV(JN  ,JM,J)=PDIV(JN  ,JM,J)-ZIN*PV(JN+1,JM,J)
        PDIV(JN+1,JM,J)=PDIV(JN+1,JM,J)+ZIN*PV(JN  ,JM,J)
      ! ELSE
        ! PVOR(JN  ,JM,J)=0._JPRB
        ! PVOR(JN+1,JM,J)=0._JPRB
        ! PDIV(JN  ,JM,J)=0._JPRB
        ! PDIV(JN+1,JM,J)=0._JPRB
      ! ENDIF
    ENDDO
  ENDDO
ENDDO
!$acc end parallel loop

IF (LHOOK) CALL DR_HOOK('EUVTVD_MOD:EUVTVD',1,ZHOOK_HANDLE)

END SUBROUTINE EUVTVD
END MODULE EUVTVD_MOD