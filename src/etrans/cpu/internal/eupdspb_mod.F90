MODULE EUPDSPB_MOD
CONTAINS
SUBROUTINE EUPDSPB(KM,KFIELD,POA,PSPEC,KFLDPTR)

!**** *EUPDSPB* - Update spectral arrays after direct Legendre transform

!     Purpose.
!     --------
!        To update spectral arrays for a fixed zonal wave-number
!         from values in POA.

!**   Interface.
!     ----------
!        CALL EUPDSPB(....)

!        Explicit arguments :  KM - zonal wavenumber
!        --------------------  KFIELD  - number of fields
!                              POA - work array
!                              PSPEC - spectral array

!        Implicit arguments :  None
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 88-02-02
!        D. Giard : 93-03-19 truncations NSMAX and NTMAX (see NOTE)
!        R. El Khatib : 94-08-02 Replace number of fields by indexes of the
!                       first and last field
!        L. Isaksen : 95-06-06 Reordering of spectral arrays
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!USE TPM_DIM
!USE TPM_FIELDS
!USE TPM_DISTR
USE TPMALD_DISTR    ,ONLY : DALD
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KM,KFIELD
REAL(KIND=JPRB)   ,INTENT(IN)  :: POA(:,:)
REAL(KIND=JPRB)   ,INTENT(OUT) :: PSPEC(:,:)
INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL :: KFLDPTR(:)

INTEGER(KIND=JPIM) :: II, INM, IR, JFLD, JN,IFLD
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.    UPDATE SPECTRAL FIELDS.
!              -----------------------
IF (LHOOK) CALL DR_HOOK('EUPDSPB_MOD:EUPDSPB',0,ZHOOK_HANDLE)
IF(PRESENT(KFLDPTR)) THEN
  DO JN=1,DALD%NCPL2M(KM),2
    INM=DALD%NESM0(KM)+(JN-1)*2
    DO JFLD=1,KFIELD
      IR= 2*JFLD-1
      II=IR+1
      IFLD = KFLDPTR(JFLD)
      PSPEC(IFLD,INM)    =POA(JN,IR)
      PSPEC(IFLD,INM+1)  =POA(JN+1,IR)
      PSPEC(IFLD,INM+2)  =POA(JN,II)
      PSPEC(IFLD,INM+3)  =POA(JN+1,II)
    ENDDO
  ENDDO
ELSE
  DO JN=1,DALD%NCPL2M(KM),2
    INM=DALD%NESM0(KM)+(JN-1)*2
! use unroll to provoke vectorization of outer loop
!cdir unroll=4
!DIR$ IVDEP
!OCL NOVREC
    DO JFLD=1,KFIELD
      IR= 2*JFLD-1
      II=IR+1
      PSPEC(JFLD,INM)    =POA(JN,IR)
      PSPEC(JFLD,INM+1)  =POA(JN+1,IR)
      PSPEC(JFLD,INM+2)  =POA(JN,II)
      PSPEC(JFLD,INM+3)  =POA(JN+1,II)
    ENDDO
  ENDDO
ENDIF
IF (LHOOK) CALL DR_HOOK('EUPDSPB_MOD:EUPDSPB',1,ZHOOK_HANDLE)

END SUBROUTINE EUPDSPB
END MODULE EUPDSPB_MOD
