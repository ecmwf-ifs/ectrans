SUBROUTINE SPEC_SETUP4PY(KRETURNCODE, KSIZEI, KSIZEJ, KPHYSICALSIZEI, KPHYSICALSIZEJ, &
                       &KTRUNCX, KTRUNCY, KNUMMAXRESOL, KLOEN, LDLAM, &
                       &KSIZEKLOEN, PDELTAX, PDELTAY, &
                       &KIDENTRESOL, LDSTOP)
! ** PURPOSE
!    Setup spectral transform for LAM and global
!
! ** DUMMY ARGUMENTS
!    KRETURNCODE: error code
!    KSIZEI, KSIZEJ: size of grid-point field (with extension zone for LAM), put max size for KSIZEI in global
!    KPHYSICALSIZEI, KPHYSICALSIZEJ: size of physical part of grid-point field for LAM (put 0 for global)
!    KTRUNCX, KTRUNCY: troncatures for LAM (only KTRUNCX is used for global, put 0 for KTRUNCY)
!    KNUMMAXRESOL: maximum number of troncatures handled
!    KLOEN: number of points on each latitude row
!    KSIZEKLOEN: size of KLOEN array
!    PDELTAX: x resolution
!    PDELTAY: y resolution
!    LDLAM: LAM (.TRUE.) or global (.FALSE.)
!    KIDENTRESOL: identification of resolution
!    LDSTOP: exception raised?
!
! ** AUTHOR
!    9 April 2014, S. Riette
!
! ** MODIFICATIONS
!    6 Jan 2016, S. Riette: PDELTAX and PDELTAY added
!    31 Jan 2019 R. El Khatib fix for single precision compilation
!
! I. Dummy arguments declaration
USE PARKIND1, ONLY : JPRB
IMPLICIT NONE
INTEGER(KIND=8), INTENT(OUT) :: KRETURNCODE
INTEGER, INTENT(IN) :: KSIZEI, KSIZEJ
INTEGER, INTENT(IN) :: KPHYSICALSIZEI, KPHYSICALSIZEJ
INTEGER, INTENT(IN) :: KTRUNCX, KTRUNCY
INTEGER, INTENT(IN) :: KNUMMAXRESOL
INTEGER, DIMENSION(KSIZEKLOEN), INTENT(IN) :: KLOEN
LOGICAL, INTENT(IN) :: LDLAM
INTEGER, INTENT(IN) :: KSIZEKLOEN
REAL(KIND=8), INTENT(IN) :: PDELTAX
REAL(KIND=8), INTENT(IN) :: PDELTAY
INTEGER, INTENT(OUT) :: KIDENTRESOL
LOGICAL, INTENT(OUT) :: LDSTOP
!
! II. Local variables declaration
INTEGER, DIMENSION(2*KSIZEJ) :: ILOEN
INTEGER :: JI
LOGICAL, SAVE :: LLFIRSTCALL=.TRUE.
LOGICAL :: LLNEWRESOL
INTEGER, SAVE :: INBRESOL=0
INTEGER(KIND=8) :: ICODEILOEN
INTEGER, SAVE :: INUMMAXRESOLSAVE=-1
INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: ITRUNCXSAVE, ITRUNCYSAVE, &
                                            IPHYSICALSIZEISAVE, &
                                            IPHYSICALSIZEJSAVE, &
                                            ISIZEISAVE, ISIZEJSAVE, &
                                            IIDENTRESOLSAVE
INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: ILOENSAVE
REAL(KIND=8), DIMENSION(:), ALLOCATABLE, SAVE :: ZDELTAXSAVE, &
                                                        ZDELTAYSAVE
REAL(KIND=8) :: ZEXWN, ZEYWN

#include "setup_trans0.h"
#include "esetup_trans.h"
#include "setup_trans.h"

KRETURNCODE=0
LDSTOP=.FALSE.
! III. Setup

! III.a Setup LAM and global spectral transform - all resolutions
! Maximum number of resolution is set now and cannot be change anymore
IF (LLFIRSTCALL) THEN
  !This code is called only once, whatever is the number of resolutions
  CALL SETUP_TRANS0(KPRINTLEV=0, LDMPOFF=.TRUE., KMAX_RESOL=KNUMMAXRESOL)
  ALLOCATE(ITRUNCXSAVE(KNUMMAXRESOL))
  ALLOCATE(ITRUNCYSAVE(KNUMMAXRESOL))
  ALLOCATE(IPHYSICALSIZEISAVE(KNUMMAXRESOL))
  ALLOCATE(IPHYSICALSIZEJSAVE(KNUMMAXRESOL))
  ALLOCATE(ISIZEJSAVE(KNUMMAXRESOL))
  ALLOCATE(ISIZEISAVE(KNUMMAXRESOL))
  ALLOCATE(ILOENSAVE(KNUMMAXRESOL))
  ALLOCATE(IIDENTRESOLSAVE(KNUMMAXRESOL))
  ALLOCATE(ZDELTAXSAVE(KNUMMAXRESOL))
  ALLOCATE(ZDELTAYSAVE(KNUMMAXRESOL))
  ITRUNCXSAVE=-1
  ITRUNCYSAVE=-1
  IPHYSICALSIZEISAVE=-1
  IPHYSICALSIZEJSAVE=-1
  ISIZEJSAVE=-1
  ISIZEISAVE=-1
  ILOENSAVE=-1
  IIDENTRESOLSAVE=-1
  ZDELTAXSAVE=-1.
  ZDELTAXSAVE=-1.
  LLFIRSTCALL=.FALSE.
  INUMMAXRESOLSAVE=KNUMMAXRESOL
ENDIF
!
! III.b Is-it a new resolution?
LLNEWRESOL=.TRUE.
IF(LDLAM) THEN
  ILOEN(:)=KSIZEI
ELSE
  ILOEN(:)=0
  ILOEN(1:MIN(SIZE(ILOEN),SIZE(KLOEN)))=KLOEN(1:MIN(SIZE(ILOEN),SIZE(KLOEN)))
ENDIF
ICODEILOEN=0
DO JI=1, SIZE(ILOEN)
  ICODEILOEN=ICODEILOEN+ILOEN(JI)*JI**4
ENDDO
DO JI=1, INBRESOL
  IF (KTRUNCX==ITRUNCXSAVE(JI) .AND. KTRUNCY==ITRUNCYSAVE(JI) .AND. &
   &KPHYSICALSIZEI==IPHYSICALSIZEISAVE(JI) .AND. &
   &KPHYSICALSIZEJ==IPHYSICALSIZEJSAVE(JI) .AND. &
   &KSIZEJ==ISIZEJSAVE(JI) .AND. KSIZEI==ISIZEISAVE(JI) .AND. &
   &ICODEILOEN==ILOENSAVE(JI) .AND. &
   &PDELTAX==ZDELTAXSAVE(JI) .AND. PDELTAY==ZDELTAYSAVE(JI)) THEN
    KIDENTRESOL=IIDENTRESOLSAVE(JI)
    LLNEWRESOL=.FALSE.
  ENDIF
ENDDO
IF(LLNEWRESOL) THEN
  INBRESOL=INBRESOL+1
  IF(INBRESOL>INUMMAXRESOLSAVE) THEN
    PRINT*, "Error in SPEC_SETUP4PY : Maximum number of resolution is exceeded."
    KRETURNCODE=-999
    LDSTOP=.TRUE.
  ENDIF
ENDIF
!
! III.c Setup LAM or global spectral transform - once by resolution
IF(LLNEWRESOL .AND. .NOT. LDSTOP) THEN
  ! The following code is exectuded once for each resolution
  ITRUNCXSAVE(INBRESOL)=KTRUNCX
  ITRUNCYSAVE(INBRESOL)=KTRUNCY
  IPHYSICALSIZEISAVE(INBRESOL)=KPHYSICALSIZEI
  IPHYSICALSIZEJSAVE(INBRESOL)=KPHYSICALSIZEJ
  ISIZEISAVE(INBRESOL)=KSIZEI
  ISIZEJSAVE(INBRESOL)=KSIZEJ
  ILOENSAVE(INBRESOL)=ICODEILOEN
  ZDELTAXSAVE(INBRESOL)=PDELTAX
  ZDELTAYSAVE(INBRESOL)=PDELTAY
  IF(LDLAM) THEN
    ZEXWN=2*3.141592653589797/(KSIZEI*PDELTAX)
    ZEYWN=2*3.141592653589797/(KSIZEJ*PDELTAY)
    CALL ESETUP_TRANS(KMSMAX=ITRUNCXSAVE(INBRESOL), KSMAX=ITRUNCYSAVE(INBRESOL), &
                     &KDGUX=IPHYSICALSIZEJSAVE(INBRESOL), &
                     &KDGL=ISIZEJSAVE(INBRESOL), KLOEN=ILOEN(:), KRESOL=IIDENTRESOLSAVE(INBRESOL), &
                     &PEXWN=REAL(ZEXWN,KIND=JPRB), PEYWN=REAL(ZEYWN,KIND=JPRB))
  ELSE
    PRINT*, "Setup spectral transform"
    CALL SETUP_TRANS(KSMAX=ITRUNCXSAVE(INBRESOL), KDGL=ISIZEJSAVE(INBRESOL), &
                    &KLOEN=ILOEN(1:ISIZEJSAVE(INBRESOL)), KRESOL=IIDENTRESOLSAVE(INBRESOL))
    PRINT*, "End Setup spectral transform"
  ENDIF
  KIDENTRESOL=IIDENTRESOLSAVE(INBRESOL)
ENDIF
END SUBROUTINE SPEC_SETUP4PY

