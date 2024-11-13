SUBROUTINE SP2GP_LAM4PY(KRETURNCODE, KSIZEI, KSIZEJ, KPHYSICALSIZEI, KPHYSICALSIZEJ, &
                         &KTRUNCX, KTRUNCY, KNUMMAXRESOL, KSIZE, LGRADIENT, LREORDER, PDELTAX, PDELTAY, &
                         &PSPEC, PGPT, PGPTM, PGPTL)
! ** PURPOSE
!    Transform spectral coefficients into grid-point values
!
! ** DUMMY ARGUMENTS
!    KRETURNCODE: error code
!    KSIZEI, KSIZEJ: size of grid-point field (with extension zone)
!    KPHYSICALSIZEI, KPHYSICALSIZEJ: size of physical part of grid-point field
!    KTRUNCX, KTRUNCY: troncatures
!    KNUMMAXRESOL: maximum number of troncatures handled
!    KSIZE: size of PSPEC
!    LREORDER: switch to reorder spectral coefficients or not
!    LGRADIENT: switch to compute or not gradient
!    PDELTAX: x resolution
!    PDELTAY: y resolution
!    PSPEC: spectral coefficient array
!    PGPT: grid-point field
!    PGPTM: N-S derivative if LGRADIENT
!    PGPTL: E-W derivative if LGRADIENT
!
! ** AUTHOR
!    9 April 2014, S. Riette
!
! ** MODIFICATIONS
!    5 Jan., S. Riette: PDELTAX, PDELTAY, LGRADIENT, PGPTM and PGPTL added
!    March, 2016, A.Mary: LREORDER
!
! I. Dummy arguments declaration
USE PARKIND1, ONLY : JPRB
IMPLICIT NONE
INTEGER(KIND=8), INTENT(OUT) :: KRETURNCODE
INTEGER(KIND=8), INTENT(IN) :: KSIZEI, KSIZEJ
INTEGER(KIND=8), INTENT(IN) :: KPHYSICALSIZEI, KPHYSICALSIZEJ
INTEGER(KIND=8), INTENT(IN) :: KTRUNCX, KTRUNCY
INTEGER(KIND=8), INTENT(IN) :: KNUMMAXRESOL
INTEGER(KIND=8), INTENT(IN) :: KSIZE
LOGICAL, INTENT(IN) :: LGRADIENT
LOGICAL, INTENT(IN) :: LREORDER
REAL(KIND=8), INTENT(IN) :: PDELTAX
REAL(KIND=8), INTENT(IN) :: PDELTAY
REAL(KIND=8), DIMENSION(KSIZE), INTENT(IN) :: PSPEC
REAL(KIND=8), DIMENSION(KSIZEI*KSIZEJ), INTENT(OUT) :: PGPT
REAL(KIND=8), DIMENSION(KSIZEI*KSIZEJ), INTENT(OUT) :: PGPTM
REAL(KIND=8), DIMENSION(KSIZEI*KSIZEJ), INTENT(OUT) :: PGPTL
!
! II. Local variables declaration
INTEGER, DIMENSION(0:KTRUNCX) :: IESM0
INTEGER :: IGPTOT, ISPEC
INTEGER, DIMENSION(0:KTRUNCY) :: ISPECINI, ISPECEND
REAL(KIND=8), DIMENSION(1, KSIZE) :: ZSPBUF
REAL(KIND=JPRB), DIMENSION(KSIZEI*KSIZEJ, 3, 1) :: ZGPBUF
INTEGER :: JI, JM, JN, IINDEX, IIDENTRESOL
LOGICAL :: LLSTOP
INTEGER :: ISIZEI, ISIZEJ, &
         & IPHYSICALSIZEI, IPHYSICALSIZEJ, &
         & ITRUNCX, ITRUNCY, &
         & INUMMAXRESOL
INTEGER, DIMENSION(1) :: ILOEN

#include "einv_trans.h"
#include "etrans_inq.h"

KRETURNCODE=0
LLSTOP=.FALSE.
ISIZEI=KSIZEI
ISIZEJ=KSIZEJ
IPHYSICALSIZEI=KPHYSICALSIZEI
IPHYSICALSIZEJ=KPHYSICALSIZEJ
ITRUNCX=KTRUNCX
ITRUNCY=KTRUNCY
INUMMAXRESOL=KNUMMAXRESOL
ILOEN(:)=0

! III. Setup
CALL SPEC_SETUP4PY(KRETURNCODE, ISIZEI, ISIZEJ, IPHYSICALSIZEI, IPHYSICALSIZEJ, &
                  &ITRUNCX, ITRUNCY, INUMMAXRESOL, ILOEN, .TRUE., 1, &
                  &PDELTAX, PDELTAY, IIDENTRESOL, LLSTOP)

! IV. Transformation

! IV.a Shape of coefficient array
!IGPTOT is the total number of points in grid-point space
!ISPEC is the number of spectral coefficients
!IESM0(m) is the index of spectral coefficient (m,0) in model
!ISPECINI(n) is the index of the first of the 4 spectral coefficient (0,n) in FA file
!ISPECEND(n) is the index of the last of the last 4 spectral coefficients (:,n) in FA file
IF (.NOT. LLSTOP) THEN
  CALL ETRANS_INQ(KRESOL=IIDENTRESOL, KGPTOT=IGPTOT, KSPEC=ISPEC, KESM0=IESM0)
  JI=1
  DO JN=0, ITRUNCY
    ISPECINI(JN)=(JI-1)*4+1
    JI=JI+COUNT(IESM0(1:ITRUNCX)-IESM0(0:ITRUNCX-1)>JN*4)
    IF (ISPEC-IESM0(ITRUNCX)>JN*4) JI=JI+1
    ISPECEND(JN)=(JI-1)*4
  ENDDO
ENDIF

! III.b Reordering
! reorder Aladin :  file ordering = coeffs per blocks of m, 4 reals per coeff
!           Aladin array ordering = coeffs per blocks of n, 4 reals per coeff
IF (LREORDER) THEN
  IF (.NOT. LLSTOP) THEN
    ZSPBUF(:,:)=0.
    JI=1
    DO JM=0,ITRUNCX+1
      DO JN=0,ITRUNCY
        IF (ISPECINI(JN)+JM*4+3<=ISPECEND(JN)) THEN
          DO IINDEX=ISPECINI(JN)+JM*4, ISPECINI(JN)+JM*4+3
            ZSPBUF(1,JI)=PSPEC(IINDEX)
            JI=JI+1
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    IF (JI/=ISPEC+1) THEN
      PRINT*, "Internal error in SP2GP_LAM4PY (spectral reordering)"
      KRETURNCODE=-999
      LLSTOP=.TRUE.
    ENDIF
  ENDIF
ELSE
  ZSPBUF(1,:) = PSPEC(:)
ENDIF

! III.c Inverse transform
IF (.NOT. LLSTOP) THEN
  IF (.NOT. LGRADIENT) THEN
      CALL EINV_TRANS(PSPSCALAR=REAL(ZSPBUF(:,:),KIND=JPRB), PGP=ZGPBUF(:,:,:), KRESOL=IIDENTRESOL)
      PGPT(:)=REAL(ZGPBUF(:,1,1),KIND=8)
  ELSE
      CALL EINV_TRANS(PSPSCALAR=REAL(ZSPBUF(:,:),KIND=JPRB), PGP=ZGPBUF(:,:,:), KRESOL=IIDENTRESOL, LDSCDERS=.TRUE.)
      PGPT(:)=REAL(ZGPBUF(:,1,1),KIND=8)
      PGPTM(:)=REAL(ZGPBUF(:,2,1),KIND=8)
      PGPTL(:)=REAL(ZGPBUF(:,3,1),KIND=8)
  ENDIF
ENDIF

END SUBROUTINE SP2GP_LAM4PY
