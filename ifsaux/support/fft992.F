C     SUBROUTINE 'FFT992' - MULTIPLE FAST REAL PERIODIC TRANSFORM
C
C     Author: Clive Temperton, January 1998
C         
C     This routine is a modernized and enhanced version of FFT991
C         - Cray directives and ancient Fortran constructs removed
C         - "vector chopping" removed
C         - WORK array is now dynamically allocated
C         - stride in WORK array is now always 1
C
C     REAL TRANSFORM OF LENGTH N PERFORMED BY REMOVING REDUNDANT
C     OPERATIONS FROM COMPLEX TRANSFORM OF LENGTH N
C
C     A IS THE ARRAY CONTAINING INPUT & OUTPUT DATA
C     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
C     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N
C     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
C         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
C     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
C     N IS THE LENGTH OF THE DATA VECTORS
C     LOT IS THE NUMBER OF DATA VECTORS
C     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
C           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
C
C     ORDERING OF COEFFICIENTS:
C         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
C         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
C
C     ORDERING OF DATA:
C         X(0),X(1),X(2),...,X(N-1), 0 , 0 ; (N+2) LOCATIONS REQUIRED
C
C     VECTORIZATION IS ACHIEVED BY DOING THE TRANSFORMS IN PARALLEL
C
C     N MUST BE COMPOSED OF FACTORS 2,3 & 5 BUT DOES NOT HAVE TO BE EVEN
C
C     DEFINITION OF TRANSFORMS:
C     -------------------------
C
C     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
C         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
C
C     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
C               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
C
      SUBROUTINE FFT992(A,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
      REAL A(*),TRIGS(N)
      INTEGER IFAX(10)
C     Dynamically allocated work array:
      REAL WORK(N*LOT+1)
      LOGICAL LIPL
*
      NFAX=IFAX(1)
      IF (ISIGN.EQ.+1) THEN
C
C     ISIGN=+1, SPECTRAL TO GRIDPOINT TRANSFORM
C     -----------------------------------------
C
        I=1
!OCL NOVREC
        DO J=1,LOT 
          A(I+INC)=0.5*A(I)
          I=I+JUMP
        ENDDO
        IF (MOD(N,2).EQ.0) THEN
          I=N*INC+1
!OCL NOVREC
          DO J=1,LOT 
            A(I)=0.5*A(I)
            I=I+JUMP
          ENDDO
        ENDIF
C
        IA=INC+1
        LA=1
        IGO=+1
C
        DO K=1,NFAX
          IFAC=IFAX(K+1)
          IERR=-1
          IF (K.EQ.NFAX.AND.NFAX.GT.2.AND.IGO.EQ.+1) THEN
            LIPL=.TRUE.
          ELSE
            LIPL=.FALSE.
          ENDIF  
          IF (INC.EQ.1.AND.JUMP.LT.(2*N).AND.
     *        K.GT.1.AND.K.LT.(NFAX-MOD(NFAX,2))) THEN
            INCA=LOT
            JUMPA=1
          ELSE
            INCA=INC
            JUMPA=JUMP
          ENDIF
          IF (IGO.EQ.+1) THEN
           CALL RPASSF(A(IA),A(IA+LA*INCA),WORK(1),WORK(IFAC*LA*LOT+1),
     *                  TRIGS,INCA,LOT,JUMPA,1,LOT,N,IFAC,LA,IERR,LIPL)
          ELSE
           CALL RPASSF(WORK(1),WORK(LA*LOT+1),A(IA),A(IA+IFAC*LA*INCA),
     *                  TRIGS,LOT,INCA,1,JUMPA,LOT,N,IFAC,LA,IERR,LIPL)
          ENDIF
          IF (IERR.NE.0) THEN
            IF (IERR.EQ.2) WRITE(6,901) IFAC
            IF (IERR.EQ.3) WRITE(6,902) IFAC
            RETURN
          ENDIF
          LA=IFAC*LA
          IGO=-IGO
          IA=1
        ENDDO
C
C     IF NECESSARY, COPY RESULTS BACK TO A
C     ------------------------------------
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
C
C     FILL IN ZEROS AT END
C     --------------------
        IX=N*INC+1
!OCL NOVREC
        DO J=1,LOT 
          A(IX)=0.0
          A(IX+INC)=0.0
          IX=IX+JUMP
        ENDDO
C
      ELSEIF (ISIGN.EQ.-1) THEN
C
C     ISIGN=-1, GRIDPOINT TO SPECTRAL TRANSFORM
C     -----------------------------------------
        IA=1
        LA=N
        IGO=+1
C
        DO K=1,NFAX
          IFAC=IFAX(NFAX+2-K)
          LA=LA/IFAC
          IERR=-1
          IF (K.EQ.1.AND.NFAX.GT.2.AND.MOD(NFAX,2).EQ.1) THEN
            LIPL=.TRUE.
          ELSE
            LIPL=.FALSE.
          ENDIF  
          IF (INC.EQ.1.AND.JUMP.LT.(2*N).AND.
     *        K.GT.(1+MOD(NFAX,2)).AND.K.LT.NFAX) THEN
            INCA=LOT
            JUMPA=1
          ELSE
            INCA=INC
            JUMPA=JUMP
          ENDIF
          IF (IGO.EQ.+1) THEN
           CALL QPASSF(A(IA),A(IA+IFAC*LA*INCA),WORK(1),WORK(LA*LOT+1),
     *                  TRIGS,INCA,LOT,JUMPA,1,LOT,N,IFAC,LA,IERR,LIPL)
          ELSE
           CALL QPASSF(WORK(1),WORK(IFAC*LA*LOT+1),A(IA),A(IA+LA*INCA),
     *                  TRIGS,LOT,INCA,1,JUMPA,LOT,N,IFAC,LA,IERR,LIPL)
          ENDIF
          IF (IERR.NE.0) THEN
            IF (IERR.EQ.2) WRITE(6,901) IFAC
            IF (IERR.EQ.3) WRITE(6,902) IFAC
            RETURN
          ENDIF
          IF (LIPL) THEN
            IA=1
          ELSE
            IGO=-IGO
            IA=INC+1
          ENDIF
        ENDDO
C
C     IF NECESSARY, COPY RESULTS BACK TO A
C     ------------------------------------
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
C
C     SHIFT A(0) & FILL IN ZERO IMAG PARTS
C     ------------------------------------
        IX=1
!OCL NOVREC
        DO J=1,LOT
          A(IX)=A(IX+INC)
          A(IX+INC)=0.0
          IX=IX+JUMP
        ENDDO
        IF (MOD(N,2).EQ.0) THEN
          IX=(N+1)*INC+1
          DO J=1,LOT
            A(IX)=0.0
            IX=IX+JUMP
          ENDDO
        ENDIF
C
      ENDIF
C
      RETURN
C
C     FORMAT STATEMENTS FOR ERROR MESSAGES:
  901 FORMAT(' FACTOR =',I3,' NOT CATERED FOR')
  902 FORMAT(' FACTOR =',I3,' ONLY CATERED FOR IF LA*IFAC=N')
C
      END
