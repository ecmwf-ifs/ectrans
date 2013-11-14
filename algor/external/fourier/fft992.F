!
!     SUBROUTINE 'FFT992' - MULTIPLE FAST REAL PERIODIC TRANSFORM
!
!     Author: Clive Temperton, January 1998
!
!     This routine is a modernized and enhanced version of FFT991
!         - Cray directives and ancient Fortran constructs removed
!         - "vector chopping" removed
!         - WORK array is now dynamically allocated
!         - stride in WORK array is now always 1
!
!     REAL TRANSFORM OF LENGTH N PERFORMED BY REMOVING REDUNDANT
!     OPERATIONS FROM COMPLEX TRANSFORM OF LENGTH N
!
!     A IS THE ARRAY CONTAINING INPUT & OUTPUT DATA
!     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
!     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N
!     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
!         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
!     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
!     N IS THE LENGTH OF THE DATA VECTORS
!     LOT IS THE NUMBER OF DATA VECTORS
!     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
!           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
!
!     ORDERING OF COEFFICIENTS:
!         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
!         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
!
!     ORDERING OF DATA:
!         X(0),X(1),X(2),...,X(N-1), 0 , 0 ; (N+2) LOCATIONS REQUIRED
!
!     VECTORIZATION IS ACHIEVED BY DOING THE TRANSFORMS IN PARALLEL
!
!     N MUST BE COMPOSED OF FACTORS 2,3 & 5 BUT DOES NOT HAVE TO BE EVEN
!
!     DEFINITION OF TRANSFORMS:
!     -------------------------
!
!     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
!         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
!
!     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
!               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
!
#ifdef MATHKEISAN
      SUBROUTINE FFT992(A,TRIGS_,IFAX_,INC,JUMP,N,LOT,ISIGN)
!AUTOPROMOTE
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      implicit none
      integer :: n
      REAL A(*),TRIGS_(N)
      INTEGER IFAX_(10)

      integer :: inc, jump, lot, isign

      real ,allocatable,dimension(:),save :: work , trigs
      integer,save :: ifax (32)


      integer, save ::  n_old=-1
      integer, save ::  lot_old=-1

!$OMP threadprivate(ifax,n_old,lot_old,trigs,work)


      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FFT992',0,ZHOOK_HANDLE)
      if (n .ne. n_old) then

         IF( ALLOCATED( work ) ) DEALLOCATE( work )
         IF( ALLOCATED( trigs ) ) DEALLOCATE( trigs )

         allocate(work(3*n*lot))
         allocate(trigs(2*n))

         call dftfax ( n, ifax, trigs )

         n_old = n
         lot_old = lot

      else

       if (lot .gt. lot_old) then

         IF( ALLOCATED( work ) ) DEALLOCATE( work )
         allocate(work(3*n*lot))
         lot_old = lot

       endif

      endif

      call dfftmlt ( a, work, trigs, ifax, inc, jump, n, lot, isign )

      IF (LHOOK) CALL DR_HOOK('FFT992',1,ZHOOK_HANDLE)
      return

      END SUBROUTINE FFT992
#else
      SUBROUTINE FFT992(A,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
!AUTOPROMOTE
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
!      
      INTEGER :: NU,N,INC,JBASE,JUMP,J,JJ,JUMPA,LOT,K,LA,NFAX
      INTEGER :: ISIGN,I,IA,IBASE,IERR,IFAC,IGO,II,INCA,IX

      REAL(KIND=8) :: A(*),TRIGS(N)
      INTEGER IFAX(10)
!     Dynamically allocated work array:
      REAL(KIND=8) :: WORK(N*LOT+1)
      LOGICAL LIPL
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FFT992',0,ZHOOK_HANDLE)
      NFAX=IFAX(1)
      IF (ISIGN.EQ.+1) THEN
!
!     ISIGN=+1, SPECTRAL TO GRIDPOINT TRANSFORM
!     -----------------------------------------
!
        I=1
!OCL NOVREC
        DO J=1,LOT
          A(I+INC)=0.5_8*A(I)
          I=I+JUMP
        ENDDO
        IF (MOD(N,2).EQ.0) THEN
          I=N*INC+1
!OCL NOVREC
          DO J=1,LOT
            A(I)=0.5_8*A(I)
            I=I+JUMP
          ENDDO
        ENDIF
!
        IA=INC+1
        LA=1
        IGO=+1
!
        DO K=1,NFAX
          IFAC=IFAX(K+1)
          IERR=-1
          IF (K.EQ.NFAX.AND.NFAX.GT.2.AND.IGO.EQ.+1) THEN
            LIPL=.TRUE.
          ELSE
            LIPL=.FALSE.
          ENDIF
          IF (INC.EQ.1.AND.JUMP.LT.(2*N).AND.                           &
     &        K.GT.1.AND.K.LT.(NFAX-MOD(NFAX,2))) THEN
            INCA=LOT
            JUMPA=1
          ELSE
            INCA=INC
            JUMPA=JUMP
          ENDIF
          IF (IGO.EQ.+1) THEN
           CALL RPASSF(A(IA),A(IA+LA*INCA),WORK(1),WORK(IFAC*LA*LOT+1), &
     &                  TRIGS,INCA,LOT,JUMPA,1,LOT,N,IFAC,LA,IERR,LIPL)
          ELSE
           CALL RPASSF(WORK(1),WORK(LA*LOT+1),A(IA),A(IA+IFAC*LA*INCA), &
     &                  TRIGS,LOT,INCA,1,JUMPA,LOT,N,IFAC,LA,IERR,LIPL)
          ENDIF
          IF (IERR.NE.0) THEN
            IF (IERR.EQ.2) WRITE(6,901) IFAC
            IF (IERR.EQ.3) WRITE(6,902) IFAC
            IF (LHOOK) CALL DR_HOOK('FFT992',1,ZHOOK_HANDLE)
            RETURN
          ENDIF
          LA=IFAC*LA
          IGO=-IGO
          IA=1
        ENDDO
!
!     IF NECESSARY, COPY RESULTS BACK TO A
!     ------------------------------------
        IF (NFAX.EQ.1) THEN
          IBASE=1
          JBASE=1
          DO JJ=1,N
            I=IBASE
            J=JBASE
            DO II=1,LOT
              A(J)=WORK(I)
              I=I+1
              J=J+JUMP
            ENDDO
            IBASE=IBASE+LOT
            JBASE=JBASE+INC
          ENDDO
        ENDIF
!
!     FILL IN ZEROS AT END
!     --------------------
        IX=N*INC+1
!OCL NOVREC
        DO J=1,LOT
          A(IX)=0.0_8
          A(IX+INC)=0.0_8
          IX=IX+JUMP
        ENDDO
!
      ELSEIF (ISIGN.EQ.-1) THEN
!
!     ISIGN=-1, GRIDPOINT TO SPECTRAL TRANSFORM
!     -----------------------------------------
        IA=1
        LA=N
        IGO=+1
!
        DO K=1,NFAX
          IFAC=IFAX(NFAX+2-K)
          LA=LA/IFAC
          IERR=-1
          IF (K.EQ.1.AND.NFAX.GT.2.AND.MOD(NFAX,2).EQ.1) THEN
            LIPL=.TRUE.
          ELSE
            LIPL=.FALSE.
          ENDIF
          IF (INC.EQ.1.AND.JUMP.LT.(2*N).AND.                           &
     &        K.GT.(1+MOD(NFAX,2)).AND.K.LT.NFAX) THEN
            INCA=LOT
            JUMPA=1
          ELSE
            INCA=INC
            JUMPA=JUMP
          ENDIF
          IF (IGO.EQ.+1) THEN
           CALL QPASSF(A(IA),A(IA+IFAC*LA*INCA),WORK(1),WORK(LA*LOT+1), &
     &                  TRIGS,INCA,LOT,JUMPA,1,LOT,N,IFAC,LA,IERR,LIPL)
          ELSE
           CALL QPASSF(WORK(1),WORK(IFAC*LA*LOT+1),A(IA),A(IA+LA*INCA), &
     &                  TRIGS,LOT,INCA,1,JUMPA,LOT,N,IFAC,LA,IERR,LIPL)
          ENDIF
          IF (IERR.NE.0) THEN
            IF (IERR.EQ.2) WRITE(6,901) IFAC
            IF (IERR.EQ.3) WRITE(6,902) IFAC
            IF (LHOOK) CALL DR_HOOK('FFT992',1,ZHOOK_HANDLE)
            RETURN
          ENDIF
          IF (LIPL) THEN
            IA=1
          ELSE
            IGO=-IGO
            IA=INC+1
          ENDIF
        ENDDO
!
!     IF NECESSARY, COPY RESULTS BACK TO A
!     ------------------------------------
        IF (NFAX.EQ.1) THEN
          IBASE=1
          JBASE=INC+1
          DO JJ=1,N
            I=IBASE
            J=JBASE
            DO II=1,LOT
              A(J)=WORK(I)
              I=I+1
              J=J+JUMP
            ENDDO
            IBASE=IBASE+LOT
            JBASE=JBASE+INC
          ENDDO
        ENDIF
!
!     SHIFT A(0) & FILL IN ZERO IMAG PARTS
!     ------------------------------------
        IX=1
!OCL NOVREC
        DO J=1,LOT
          A(IX)=A(IX+INC)
          A(IX+INC)=0.0_8
          IX=IX+JUMP
        ENDDO
        IF (MOD(N,2).EQ.0) THEN
          IX=(N+1)*INC+1
          DO J=1,LOT
            A(IX)=0.0_8
            IX=IX+JUMP
          ENDDO
        ENDIF
!
      ENDIF
!
!     FORMAT STATEMENTS FOR ERROR MESSAGES:
  901 FORMAT(' FACTOR =',I3,' NOT CATERED FOR')
  902 FORMAT(' FACTOR =',I3,' ONLY CATERED FOR IF LA*IFAC=N')
!
      IF (LHOOK) CALL DR_HOOK('FFT992',1,ZHOOK_HANDLE)
      ENDSUBROUTINE FFT992
#endif
