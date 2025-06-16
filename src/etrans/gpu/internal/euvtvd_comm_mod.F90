MODULE EUVTVD_COMM_MOD
CONTAINS
SUBROUTINE EUVTVD_COMM(KM,KMLOC,KFIELD,KFLDPTR,PU,PV,PSPMEANU,PSPMEANV)

!**** *EUVTVD* - Compute vor/div from u and v in spectral space

!     Purpose.
!     --------
!        To compute vorticity and divergence from u and v in spectral
!       space. Input u and v from KM to NTMAX+1, output vorticity and
!       divergence from KM to NTMAX - communication part.

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
!        N. Lopes & R. El Khatib 15-Jun-2012 Scalability enhancement
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DIM
USE TPM_FIELDS
USE TPM_DISTR
USE TPMALD_GEO
USE TPMALD_DISTR
USE MPL_MODULE
USE SET2PE_MOD
USE ABORT_TRANS_MOD
IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KFIELD, KM, KMLOC
REAL(KIND=JPRB), INTENT(INOUT) :: PU  (:,:,:),PV  (:,:,:)

INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN)  :: KFLDPTR(:)
REAL(KIND=JPRB),    OPTIONAL, INTENT(OUT) :: PSPMEANU(:),PSPMEANV(:)

INTEGER(KIND=JPIM) :: II, IN, IR, J, JN

INTEGER(KIND=JPIM) :: ISENDREQ(NPRTRW)

REAL(KIND=JPRB) :: ZSPU(2*KFIELD)
REAL(KIND=JPRB) :: ZIN
INTEGER(KIND=JPIM) :: JA,ITAG,ILEN,IFLD,ISND
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EUVTVD_COMM_MOD:EUVTVD_COMM',0,ZHOOK_HANDLE)

!*       1.    COMPUTE U V FROM VORTICITY AND DIVERGENCE.
!              ------------------------------------------

IF (KM == 0) THEN

!$acc data present(PU,PV)
!$acc data copyout (PSPMEANU, PSPMEANV) copyin(KMLOC)
!$acc data copyin (KFLDPTR) if(present (KFLDPTR))

  IF (PRESENT(KFLDPTR)) THEN
!$acc parallel loop private(ir,ifld)
    DO J = 1, KFIELD
      IR = 2*J-1
      IFLD=KFLDPTR(J)
      PSPMEANU(IFLD)=PU(1,KMLOC,IR)
      PSPMEANV(IFLD)=PV(1,KMLOC,IR)
    ENDDO
!$acc end parallel loop 
  ELSE
!$acc parallel loop private(j,ir)
    DO J = 1, KFIELD
      IR = 2*J-1
      PSPMEANU(J)=PU(1,KMLOC,IR)
      PSPMEANV(J)=PV(1,KMLOC,IR)
    ENDDO
!$acc end parallel loop 
  ENDIF

!$acc end data
!$acc end data
!$acc end data
ENDIF

IF (NPRTRW > 1 .AND. KFIELD > 0) THEN
  IF (KM == 0) THEN
    IF (PRESENT(KFLDPTR)) THEN
      DO J=1,KFIELD
        IFLD=KFLDPTR(J)
        ZSPU(J)=PSPMEANU(IFLD)
        ZSPU(KFIELD+J)=PSPMEANV(IFLD)
      ENDDO 
    ELSE
      DO J=1,KFIELD
        ZSPU(J)=PSPMEANU(J)
        ZSPU(KFIELD+J)=PSPMEANV(J)
      ENDDO
    ENDIF
    DO JA=1,NPRTRW
      IF (JA /= MYSETW) THEN
        CALL SET2PE(ISND,0,0,JA,MYSETV)
        ISND=NPRCIDS(ISND)
        ITAG=300000+KFIELD*NPROC+ISND
        CALL MPL_SEND(ZSPU(1:2*KFIELD),KDEST=ISND,KTAG=ITAG, &
         &   KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(JA),CDSTRING='EUVTVD_COMM:')
      ENDIF
    ENDDO
  ELSE
    IF (KMLOC == 1) THEN
      IF (D%NPROCM(0) /= MYSETW) THEN
        CALL SET2PE(ISND,0,0,D%NPROCM(0),MYSETV)
        ISND=NPRCIDS(ISND)          
        ITAG=300000+KFIELD*NPROC+MYPROC

        CALL MPL_RECV(ZSPU(1:2*KFIELD),KSOURCE=ISND,KTAG=ITAG,KOUNT=ILEN, &
         &   CDSTRING='EUVTVD_COMM:')
        IF (ILEN /= 2*KFIELD) CALL ABORT_TRANS('EUVTVD_COMM: RECV INVALID RECEIVE MESSAGE LENGHT')
        IF (PRESENT(KFLDPTR)) THEN 
          DO J=1,KFIELD
            IFLD=KFLDPTR(J)
            PSPMEANU(IFLD)=ZSPU(J)
            PSPMEANV(IFLD)=ZSPU(KFIELD+J)
          ENDDO
        ELSE
          DO J=1,KFIELD
            PSPMEANU(J)=ZSPU(J)
            PSPMEANV(J)=ZSPU(KFIELD+J)
          ENDDO
        ENDIF
      ENDIF
    ENDIF
  ENDIF
ENDIF

IF (LHOOK) CALL DR_HOOK('EUVTVD_COMM_MOD:EUVTVD_COMM',1,ZHOOK_HANDLE)

END SUBROUTINE EUVTVD_COMM
END MODULE EUVTVD_COMM_MOD