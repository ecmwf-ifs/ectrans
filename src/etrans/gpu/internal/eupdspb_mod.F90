MODULE EUPDSPB_MOD
CONTAINS
SUBROUTINE EUPDSPB(KFIELD,POA,PSPEC,KFLDPTR)

!**** *EUPDSPB* - Update spectral arrays after direct Legendre transform

!     Purpose.
!     --------
!        To update spectral arrays for a fixed zonal wave-number
!         from values in POA.

!**   Interface.
!     ----------
!        CALL EUPDSPB(....)

!        Explicit arguments :  
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

USE TPMALD_DISTR    ,ONLY : DALD, DALD_NESM0, DALD_NCPL2M
USE TPM_DISTR       ,ONLY : D, D_MYMS, D_NUMP
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELD
REAL(KIND=JPRB)   ,INTENT(IN)  :: POA(:,:,:)
REAL(KIND=JPRB)   ,INTENT(OUT) :: PSPEC(:,:)
INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL :: KFLDPTR(:)

INTEGER(KIND=JPIM) :: II, INM, IR, JFLD, JN,IFLD, JM, IM
INTEGER(KIND=JPIM) :: JFLDPTR(KFIELD)
INTEGER(KINd=JPIM) :: MAX_NCPL2M
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


!     ------------------------------------------------------------------

!*       1.    UPDATE SPECTRAL FIELDS.
!              -----------------------
IF (LHOOK) CALL DR_HOOK('EUPDSPB_MOD:EUPDSPB',0,ZHOOK_HANDLE)

!$ACC data present (POA, PSPEC)

IF(PRESENT(KFLDPTR)) THEN
  JFLDPTR=KFLDPTR
ELSE
  DO JFLD=1,KFIELD
    JFLDPTR(JFLD)=JFLD
  ENDDO
ENDIF  

MAX_NCPL2M = MAXVAL (DALD_NCPL2M)

!$ACC parallel loop collapse(3) &
!$acc& copyin(MAX_NCPL2M,KFIELD,JFLDPTR) &
!$acc& present(D_NUMP,D_MYMS,DALD_NESM0,DALD_NCPL2M) &
!$acc& private(JM,JN,JFLD,IM,INM,IR,II )
DO JN=1,MAX_NCPL2M,2
  DO JM = 1, D_NUMP
     DO JFLD=1,KFIELD
      IM = D_MYMS(JM)
      INM=DALD_NESM0(IM)+(JN-1)*2
      if ( JN .LE. DALD_NCPL2M(IM) ) then
        IR= 2*JFLD-1
        II=IR+1
        PSPEC(JFLDPTR(JFLD),INM)    =POA(JN  ,JM,IR)
        PSPEC(JFLDPTR(JFLD),INM+1)  =POA(JN+1,JM,IR)
        PSPEC(JFLDPTR(JFLD),INM+2)  =POA(JN  ,JM,II)
        PSPEC(JFLDPTR(JFLD),INM+3)  =POA(JN+1,JM,II)
      endif
     ENDDO
   ENDDO
 
 ENDDO
 
!$ACC end data

IF (LHOOK) CALL DR_HOOK('EUPDSPB_MOD:EUPDSPB',1,ZHOOK_HANDLE)

END SUBROUTINE EUPDSPB
END MODULE EUPDSPB_MOD