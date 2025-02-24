! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE EUVTVD_COMM_MOD
CONTAINS
SUBROUTINE EUVTVD_COMM(KFIELD,PSPMEANU,PSPMEANV,KFLDPTR)

!**** *EUVTVD_COMM* - Communicate mean wind

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        CALL EUVTVD_COMM(KFIELD,PSPMEANU,PSPMEANV,KFLDPTR)

!        Explicit arguments :  
!        --------------------  KFIELD - number of fields (levels)
!                              KFLDPTR - fields pointers

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
!        R. El Khatib 12-Jan-2020 Fix missing finalization of communications
!        R. El Khatib 02-Jun-2022 Optimization/Cleaning
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

INTEGER(KIND=JPIM), INTENT(IN) :: KFIELD
REAL(KIND=JPRB),    INTENT(INOUT) :: PSPMEANU(KFIELD)
REAL(KIND=JPRB),    INTENT(INOUT) :: PSPMEANV(KFIELD)
INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN)  :: KFLDPTR(KFIELD)

INTEGER(KIND=JPIM) :: J, JA,ITAG,ILEN,IFLD,ISND, IM, JM

INTEGER(KIND=JPIM) :: ISENDREQ(NPRTRW)

REAL(KIND=JPRB) :: ZSPU(2*KFIELD)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EUVTVD_COMM_MOD:EUVTVD_COMM',0,ZHOOK_HANDLE)

!*       1.    COMMUNICATE MEAN WIND
!              ---------------------

IF (D%NPROCM(0) == MYSETW) THEN
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
      ITAG=1
      CALL MPL_SEND(ZSPU(1:2*KFIELD),KDEST=ISND,KTAG=ITAG, &
       &   KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=ISENDREQ(JA),CDSTRING='EUVTVD_COMM:')
    ENDIF
  ENDDO
  DO JA=1,NPRTRW
    IF (JA /= MYSETW) THEN
      CALL MPL_WAIT(KREQUEST=ISENDREQ(JA),CDSTRING='EUVTVD_COMM:')
    ENDIF
  ENDDO
ELSE
  CALL SET2PE(ISND,0,0,D%NPROCM(0),MYSETV)
  ITAG=1
  CALL MPL_RECV(ZSPU(1:2*KFIELD),KSOURCE=NPRCIDS(ISND),KTAG=ITAG,KOUNT=ILEN, CDSTRING='EUVTVD_COMM:')
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

IF (LHOOK) CALL DR_HOOK('EUVTVD_COMM_MOD:EUVTVD_COMM',1,ZHOOK_HANDLE)

END SUBROUTINE EUVTVD_COMM
END MODULE EUVTVD_COMM_MOD
