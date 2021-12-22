      SUBROUTINE SGTSL(KN,PC,PD,PE,PB,KINFO)
!AUTOPROMOTE
      USE PARKIND1, ONLY : JPIM, JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
      IMPLICIT NONE
      INTEGER :: IK
      INTEGER(KIND=JPIM) :: KN
      INTEGER(KIND=JPIM) :: KINFO
      REAL(KIND=JPRB) :: PC(KN)
      REAL(KIND=JPRB) :: PD(KN)
      REAL(KIND=JPRB) :: PE(KN)
      REAL(KIND=JPRB) :: PB(KN)
!
!     dgtsl given a general tridiagonal matrix and a right hand
!     side will find the solution.
!
!     on entry
!
!        kn      integer
!                is the order of the tridiagonal matrix.
!
!        pc      real(kn)
!                is the subdiagonal of the tridiagonal matrix.
!                pc(2) through pc(kn) should contain the subdiagonal.
!                on output pc is destroyed.
!
!        pd      real(kn)
!                is the diagonal of the tridiagonal matrix.
!                on output pd is destroyed.
!
!        pe      real(kn)
!                is the superdiagonal of the tridiagonal matrix.
!                pe(1) through pe(kn-1) should contain the superdiagonal.
!                on output pe is destroyed.
!
!        pb      real(kn)
!                is the right hand side vector.
!
!     on return
!
!        pb       is the solution vector.
!
!        kinfo    integer
!                = 0 normal value.
!                = k if the k-th element of the diagonal becomes
!                    exactly zero.  the subroutine returns when
!                    this is detected.
!
!     linpack. this version dated 08/14/78 .
!     jack dongarra, argonne national laboratory.
!     95-01-11 lars isaksen , ECMWF    doctorized
!
!     no externals
!     fortran abs
!
!     internal variables
!
      INTEGER(KIND=JPIM) :: JK,JB,IKP1,INM1,INM2
      REAL(KIND=JPRB) :: ZT
!     begin block permitting ...exits to 100
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('SGTSL',0,ZHOOK_HANDLE)
      KINFO = 0
      PC(1) = PD(1)
      INM1 = KN - 1
      IF (INM1 .LT. 1) GO TO 40
         PD(1) = PE(1)
         PE(1) = 0.0_JPRB
         PE(KN) = 0.0_JPRB
!
         DO 30 JK = 1, INM1
         IKP1 = JK + 1
!
!        find the largest of the two rows
!
         IF (ABS(PC(IKP1)) .LT. ABS(PC(JK))) GO TO 10
!
!           interchange row
!
            ZT = PC(IKP1)
            PC(IKP1) = PC(JK)
            PC(JK) = ZT
            ZT = PD(IKP1)
            PD(IKP1) = PD(JK)
            PD(JK) = ZT
            ZT = PE(IKP1)
            PE(IKP1) = PE(JK)
            PE(JK) = ZT
            ZT = PB(IKP1)
            PB(IKP1) = PB(JK)
            PB(JK) = ZT
   10    CONTINUE
!
!           zero elements
!
         IF (PC(JK) .NE. 0.0_JPRB) GO TO 20
            KINFO = JK
!     ............exit
            GO TO 100
   20    CONTINUE
         ZT = -PC(IKP1)/PC(JK)
         PC(IKP1) = PD(IKP1) + ZT*PD(JK)
         PD(IKP1) = PE(IKP1) + ZT*PE(JK)
         PE(IKP1) = 0.
         PB(IKP1) = PB(IKP1) + ZT*PB(JK)
   30    CONTINUE
   40    CONTINUE
         IF (PC(KN) .NE. 0.0_JPRB) GO TO 50
            KINFO = KN
            GO TO 90
   50    CONTINUE
!
!        back solve
!
         INM2 = KN - 2
         PB(KN) = PB(KN)/PC(KN)
         IF (KN .EQ. 1) GO TO 80
            PB(INM1) = (PB(INM1) - PD(INM1)*PB(KN))/PC(INM1)
            IF (INM2 .LT. 1) GO TO 70
               DO 60 JB = 1, INM2
               IK = INM2 - JB + 1
               PB(IK) = (PB(IK) - PD(IK)*PB(IK+1) - PE(IK)*PB(IK+2))    &
     &                  /PC(IK)
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
!
      IF (LHOOK) CALL DR_HOOK('SGTSL',1,ZHOOK_HANDLE)
      ENDSUBROUTINE SGTSL
