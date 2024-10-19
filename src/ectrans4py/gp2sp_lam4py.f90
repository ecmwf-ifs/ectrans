SUBROUTINE GP2SP_LAM4PY(KRETURNCODE, KSIZE, KSIZEI, KSIZEJ, KPHYSICALSIZEI, KPHYSICALSIZEJ, &
                         &KTRUNCX, KTRUNCY, KNUMMAXRESOL, PDELTAX, PDELTAY, LREORDER, PGPT, PSPEC)
! ** PURPOSE
!    Transform grid point values into spectral coefficients
!
! ** DUMMY ARGUMENTS
!    KRETURNCODE: error code
!    KSIZE: size of spectral field
!    KSIZEI, KSIZEJ: size of grid-point field (with extension zone)
!    KPHYSICALSIZEI, KPHYSICALSIZEJ: size of physical part of grid-point field
!    KTRUNCX, KTRUNCY: troncatures
!    KNUMMAXRESOL: maximum number of troncatures handled
!    PDELTAX: x resolution
!    PDELTAY: y resolution
!    LREORDER: switch to reorder spectral coefficients or not
!    PGPT: grid-point field
!    PSPEC: spectral coefficient array
!
! ** AUTHOR
!    9 April 2014, S. Riette
!
! ** MODIFICATIONS
!    6 Jan., S. Riette: PDELTAX and PDELTAY added
!    March, 2016, A.Mary: LREORDER
!
! I. Dummy arguments declaration
USE PARKIND1, ONLY : JPRB
IMPLICIT NONE
INTEGER(KIND=8), INTENT(OUT) :: KRETURNCODE
INTEGER(KIND=8), INTENT(IN) :: KSIZE, KSIZEI, KSIZEJ
INTEGER(KIND=8), INTENT(IN) :: KPHYSICALSIZEI, KPHYSICALSIZEJ
INTEGER(KIND=8), INTENT(IN) :: KTRUNCX, KTRUNCY
INTEGER(KIND=8), INTENT(IN) :: KNUMMAXRESOL
REAL(KIND=8), INTENT(IN) :: PDELTAX
REAL(KIND=8), INTENT(IN) :: PDELTAY
LOGICAL, INTENT(IN) :: LREORDER
REAL(KIND=8), DIMENSION(KSIZEI*KSIZEJ), INTENT(IN) :: PGPT
REAL(KIND=8), DIMENSION(KSIZE), INTENT(OUT) :: PSPEC
!
! II. Local variables declaration
INTEGER, DIMENSION(0:KTRUNCX) :: IESM0
INTEGER :: IGPTOT, ISPEC
INTEGER, DIMENSION(0:KTRUNCY) :: ISPECINI, ISPECEND
REAL(KIND=JPRB), DIMENSION(1, KSIZEI*KSIZEJ) :: ZSPBUF !size over-evaluated
REAL(KIND=JPRB), DIMENSION(KSIZEI*KSIZEJ, 1, 1) :: ZGPBUF
INTEGER :: JI, JM, JN, IIDENTRESOL
LOGICAL :: LLSTOP
INTEGER :: ISIZEI, ISIZEJ, &
         & IPHYSICALSIZEI, IPHYSICALSIZEJ, &
         & ITRUNCX, ITRUNCY, &
         & INUMMAXRESOL
INTEGER, DIMENSION(1) :: ILOEN

#include "edir_trans.h"
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

! III.b transform
IF (.NOT. LLSTOP) THEN
  ZGPBUF(:,1,1)=REAL(PGPT(:),KIND=JPRB)
  CALL EDIR_TRANS(PSPSCALAR=ZSPBUF(:,:), PGP=ZGPBUF(:,:,:), KRESOL=IIDENTRESOL)
ENDIF

! III.c Reordering
! reorder Aladin :  file ordering = coeffs per blocks of m, 4 reals per coeff
!           Aladin array ordering = coeffs per blocks of n, 4 reals per coeff
IF (LREORDER) THEN
  IF (.NOT. LLSTOP) THEN
    JI=1
    PSPEC(:)=0.
    DO JM=0,ITRUNCX*4+4,4
      DO JN=0,ITRUNCY
        IF (ISPECINI(JN)+JM+3<=ISPECEND(JN)) THEN
          PSPEC(ISPECINI(JN)+JM:ISPECINI(JN)+JM+3) = REAL(ZSPBUF(1,JI:JI+3),KIND=8)
          JI=JI+4
        ENDIF
      ENDDO
    ENDDO
    IF(JI/=ISPEC+1) THEN
      PRINT*, "Internal error in GP2SP_LAM4PY (spectral reordering)"
      KRETURNCODE=-999
    ENDIF
  ENDIF
ELSE
  PSPEC(1:KSIZE) = REAL(ZSPBUF(1,1:KSIZE),KIND=8)
ENDIF

END SUBROUTINE GP2SP_LAM4PY
