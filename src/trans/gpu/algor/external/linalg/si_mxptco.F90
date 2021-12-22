SUBROUTINE SI_MXPTCO(KM,KSMAX,KFLEV,KFLSUR,PF,PALPHA,PDENIM,&
 & PEPSI,PIN,POU)  

!**** *SI_MXPTCO* - Complex multiplication by a penta-diagonal matrix
!                   the only non-zero diagonals of which are:
!                   - the main diagonal.
!                   - the lower second diagonal.
!                   - the upper second diagonal.

!     Purpose.
!     --------
!       This routine is called in the SI scheme (NH model)
!       Y(n) = K0(n) * X(n) + KM2(n) * X(n-2) + KP2(n) * X(n+2)

!       1/ The complex (diagonal) coefficients K0, KM2, KP2
!          are computed from input values.
!       2/ The matrix multiplication itself is done.


!**   Interface.
!     ----------
!        *CALL* *SI_MXPTCO(...)

!        Explicit arguments :
!        --------------------
!         KM          - Zonal wavenumber                            (input)
!         KSMAX       - Truncation limit                            (input)
!         KFLEV       - Number of levels                            (input)
!         KFLSUR      - Surdimension corresponding to KFLEV         (input)
!         PF          - Help constant (2*OMEGA*BETA*DT)             (input)
!         PALPHA      - Help array                                  (input)
!         PDENIM      - Help array                                  (input)
!         PEPSI       - Help array                                  (input)
!         PIN         - known vector.                               (input)
!         POU         - unknown vector.                             (output)

!        Implicit arguments :
!        --------------------
!        none.

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        None.

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        K. YESSAD (CNRM/GMAP)

!     Modifications.
!     --------------
!        Oct. 2012, G. Mozdzynski: OMP optimisations and remove LDVECHOR 
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

!      ----------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KM 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSMAX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLSUR 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PF
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALPHA(KM:KSMAX+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDENIM(KM:KSMAX+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEPSI(KM:KSMAX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIN(KFLSUR,2,KM:KSMAX) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: POU(KFLSUR,2,KM:KSMAX) 

!      ----------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JN

REAL(KIND=JPRB) :: ZN, ZNDEN2, ZNM1, ZNP1, ZNP2, ZF2

REAL(KIND=JPRB) :: ZALPHA(KM-1:KSMAX+1)
REAL(KIND=JPRB) :: ZDENIM(KM-1:KSMAX+1)
REAL(KIND=JPRB) :: ZEPSI(KM:KSMAX+1)

REAL(KIND=JPRB) :: ZREK0(KM:KSMAX)
REAL(KIND=JPRB) :: ZIMK0(KM:KSMAX)
REAL(KIND=JPRB) :: ZREKM2(KM:KSMAX)
REAL(KIND=JPRB) :: ZIMKM2(KM:KSMAX)
REAL(KIND=JPRB) :: ZREKP2(KM:KSMAX)
REAL(KIND=JPRB) :: ZIMKP2(KM:KSMAX)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!      ----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SI_MXPTCO',0,ZHOOK_HANDLE)

!      ----------------------------------------------------------------

!*       1.    COMPUTE COEFFICIENTS K0, KM2, KP2.
!              ----------------------------------

ZF2=PF*PF

ZALPHA(KM-1)=0.0_JPRB
ZALPHA(KM:KSMAX+1)=PALPHA(KM:KSMAX+1)

ZDENIM(KM-1)=0.0_JPRB
ZDENIM(KM:KSMAX+1)=PDENIM(KM:KSMAX+1)

ZEPSI(KM:KSMAX)=PEPSI(KM:KSMAX)
ZEPSI(KSMAX+1)=0.0_JPRB

ZREK0(:)=0.0_JPRB
ZIMK0(:)=0.0_JPRB
ZREKM2(:)=0.0_JPRB
ZIMKM2(:)=0.0_JPRB
ZREKP2(:)=0.0_JPRB
ZIMKP2(:)=0.0_JPRB

! * Real part and imaginary part of K0:

DO JN=KM,KSMAX
  ZN=REAL(JN,JPRB)
  ZNDEN2=MAX(1._JPRB,ZN*ZN)
  ZNM1=REAL(JN-1,JPRB)
  ZNP1=REAL(JN+1,JPRB)
  ZNP2=REAL(JN+2,JPRB)
  ZREK0(JN)=ZF2* &
   & (ZDENIM(JN-1)*ZEPSI(JN)*ZEPSI(JN)*ZNM1*ZNP1/ZNDEN2 &
   & + ZDENIM(JN+1)*ZEPSI(JN+1)*ZEPSI(JN+1)*ZN*ZNP2/(ZNP1*ZNP1) )
  ZIMK0(JN)= - ZALPHA(JN) + ZF2* &
   & (ZALPHA(JN-1)*ZDENIM(JN-1)*ZEPSI(JN)*ZEPSI(JN)*ZNM1*ZNP1/ZNDEN2 &
   & + ZALPHA(JN+1)*ZDENIM(JN+1)*ZEPSI(JN+1)*ZEPSI(JN+1)*ZN*ZNP2/(ZNP1*ZNP1) )
ENDDO

! * Real part and imaginary part of KM2:
DO JN=KM+2,KSMAX
  ZNM1=REAL(JN-1,JPRB)
  ZNP1=REAL(JN+1,JPRB)
  ZREKM2(JN)=ZF2*ZDENIM(JN-1)*ZEPSI(JN)*ZEPSI(JN-1)*ZNP1/ZNM1
  ZIMKM2(JN)=ZALPHA(JN-1)*ZREKM2(JN)
ENDDO

! * Real part and imaginary part of KP2:
DO JN=KM,KSMAX-2
  ZN=REAL(JN,JPRB)
  ZNP2=REAL(JN+2,JPRB)
  ZREKP2(JN)=ZF2*ZDENIM(JN+1)*ZEPSI(JN+1)*ZEPSI(JN+2)*ZN/ZNP2
  ZIMKP2(JN)=ZALPHA(JN+1)*ZREKP2(JN)
ENDDO

!       ----------------------------------------------------------------

!*       2.    MATRIX MULTIPLICATION
!              ---------------------

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JN,JLEV)
DO JN=KM,KSMAX
  DO JLEV=1,KFLSUR
    POU(JLEV,1:2,JN)=0.0_JPRB
  ENDDO
ENDDO
!$OMP END PARALLEL DO

! Vectorisation on JN:

! * Add main diagonal terms:
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JN,JLEV)
DO JN=KM,KSMAX
  DO JLEV=1,KFLEV
    POU(JLEV,1,JN)=POU(JLEV,1,JN) &
     & +ZREK0(JN)*PIN(JLEV,1,JN)-ZIMK0(JN)*PIN(JLEV,2,JN)
    POU(JLEV,2,JN)=POU(JLEV,2,JN) &
     & +ZIMK0(JN)*PIN(JLEV,1,JN)+ZREK0(JN)*PIN(JLEV,2,JN)
  ENDDO
ENDDO
!$OMP END PARALLEL DO

! * Add lower non-zero diagonal terms:
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JN,JLEV)
DO JN=KM+2,KSMAX
  DO JLEV=1,KFLEV
    POU(JLEV,1,JN)=POU(JLEV,1,JN) &
     & +ZREKM2(JN)*PIN(JLEV,1,JN-2)-ZIMKM2(JN)*PIN(JLEV,2,JN-2)
    POU(JLEV,2,JN)=POU(JLEV,2,JN) &
     & +ZIMKM2(JN)*PIN(JLEV,1,JN-2)+ZREKM2(JN)*PIN(JLEV,2,JN-2)
  ENDDO
ENDDO
!$OMP END PARALLEL DO

  ! * Add upper non-zero diagonal terms:
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JN,JLEV)
DO JN=KM,KSMAX-2
  DO JLEV=1,KFLEV
    POU(JLEV,1,JN)=POU(JLEV,1,JN) &
     & +ZREKP2(JN)*PIN(JLEV,1,JN+2)-ZIMKP2(JN)*PIN(JLEV,2,JN+2)
    POU(JLEV,2,JN)=POU(JLEV,2,JN) &
     & +ZIMKP2(JN)*PIN(JLEV,1,JN+2)+ZREKP2(JN)*PIN(JLEV,2,JN+2)
  ENDDO
ENDDO
!$OMP END PARALLEL DO

! * Reset to zero the imaginary part of POU for KM=0:
IF (KM == 0) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JN,JLEV)
  DO JN=0,KSMAX
    DO JLEV=1,KFLEV
      POU(JLEV,2,JN)=0.0_JPRB
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
ENDIF

!       ----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SI_MXPTCO',1,ZHOOK_HANDLE)
END SUBROUTINE SI_MXPTCO

