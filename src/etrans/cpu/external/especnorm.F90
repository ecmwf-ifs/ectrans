! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


SUBROUTINE ESPECNORM(PSPEC,KVSET,KMASTER,KRESOL,PMET,PNORM)

!**** *ESPECNORM* - Compute global spectral norms

!     Purpose.
!     --------
!        Interface routine for computing spectral norms

!**   Interface.
!     ----------
!     CALL ESPECNORM(...)

!     Explicit arguments : All arguments optional
!     --------------------
!     PSPEC(:,:)  - Spectral array
!     KVSET(:)    - "B-Set" for each field
!     KMASTER     - processor to recieve norms
!     KRESOL      - resolution tag  which is required ,default is the
!                   first defined resulution (input)
!     PMET(:)     - metric
!     PNORM(:)    - Norms (output for processor KMASTER)

!     Method.
!     -------

!     Externals.  ESET_RESOL - set resolution
!     ----------  ESPNORM_CTL - control routine

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!ifndef INTERFACE

USE TPM_GEN         ,ONLY : NERR
!USE TPM_DIM
USE TPM_DISTR       ,ONLY : D, NPRTRV, MYSETV, MYPROC

USE ESET_RESOL_MOD  ,ONLY : ESET_RESOL
USE ESPNORM_CTL_MOD ,ONLY : ESPNORM_CTL
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS

!endif INTERFACE

IMPLICIT NONE

! Declaration of arguments

REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)    :: PSPEC(:,:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)    :: KVSET(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)    :: KMASTER
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)    :: KRESOL
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)    :: PMET(:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT)   :: PNORM(:)
!ifndef INTERFACE

INTEGER(KIND=JPIM) :: IMASTER,IFLD,IFLD_G,J
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

! Set current resolution
IF (LHOOK) CALL DR_HOOK('ESPECNORM',0,ZHOOK_HANDLE)
CALL ESET_RESOL(KRESOL)

! Set defaults
IMASTER = 1
IFLD    = 0

IF(PRESENT(KMASTER)) THEN
  IMASTER = KMASTER
ENDIF

IF(PRESENT(KVSET)) THEN
  IFLD_G = UBOUND(KVSET,1)
  DO J=1,IFLD_G
    IF(KVSET(J) > NPRTRV) THEN
      WRITE(NERR,*) 'ESPECNORM:KVSET(J) > NPRTRV ',J,KVSET(J),NPRTRV
      CALL ABORT_TRANS('ESPECNORM:KVSET TOO LONG OR CONTAINS VALUES OUTSIDE RANGE')
    ENDIF
    IF(KVSET(J) == MYSETV) THEN
      IFLD = IFLD+1
    ENDIF
  ENDDO
ELSE
  IF(PRESENT(PSPEC)) THEN
    IFLD = UBOUND(PSPEC,1)
  ENDIF
  IFLD_G = IFLD
ENDIF

IF(NPRTRV >1) THEN
  IF(IFLD > 0 .AND. .NOT. PRESENT(KVSET)) THEN
    WRITE(NERR,*)'NPRTRV >1 AND IFLD > 0 AND NOT PRESENT(KVSET)',&
     & NPRTRV,IFLD
    CALL ABORT_TRANS('ESPECNORM: SPECIFY VERTICAL SPECTRAL DISTRIBUTION!')
  ENDIF
ENDIF
IF(MYPROC == IMASTER) THEN
  IF(.NOT. PRESENT(PNORM)) THEN
    CALL ABORT_TRANS('ESPECNORM: PNORM NOT PRESENT')
  ENDIF
  IF(UBOUND(PNORM,1) < IFLD_G) THEN
    CALL ABORT_TRANS('ESPECNORM: PNORM TOO SMALL')
  ENDIF
ENDIF
IF(IFLD > 0 ) THEN
  IF(.NOT. PRESENT(PSPEC)) THEN
    CALL ABORT_TRANS('ESPECNORM: PSPEC NOT PRESENT')
  ENDIF
  IF(UBOUND(PSPEC,1) < IFLD) THEN
    CALL ABORT_TRANS('ESPECNORM: FIRST DIMENSION OF PSPEC TOO SMALL')
  ENDIF
  IF(UBOUND(PSPEC,2) < D%NSPEC2) THEN
    CALL ABORT_TRANS('ESPECNORM: FIRST DIMENSION OF PSPEC TOO SMALL')
  ENDIF
ENDIF

CALL ESPNORM_CTL(PSPEC,IFLD,IFLD_G,KVSET,IMASTER,PMET,PNORM)
IF (LHOOK) CALL DR_HOOK('ESPECNORM',1,ZHOOK_HANDLE)

!endif INTERFACE

!     ------------------------------------------------------------------

END SUBROUTINE ESPECNORM
