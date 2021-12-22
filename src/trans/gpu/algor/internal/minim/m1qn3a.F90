      SUBROUTINE M1QN3A (YDGEOMETRY,YDFIELDS,YDMTRAJ,SIMUL,PROSCA,      &
     &                   CTONB,CTCAB,N,X,F,G,DXMIN,DF1,                 &
     &                   EPSG,IMPRES,IO,MODE,NITER,NSIM,INMEMO,M,JMIN,  &
     &                   JMAX,D,GG,DIAG,AUX,ALPHA,YBAR,SBAR,YDVARBC)
!AUTOPROMOTE
      USE PARKIND1, ONLY : JPIM, JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!----
!
!     Code d'optimisation proprement dit.
!
!----
      USE GEOMETRY_MOD , ONLY : GEOMETRY
      USE FIELDS_MOD , ONLY : FIELDS
      USE MTRAJ_MOD  , ONLY : MTRAJ
      USE VARBC_CLASS,ONLY: CLASS_VARBC
      USE CONTROL_VECTORS_MOD
!
      IMPLICIT NONE
!
!         arguments
!
      LOGICAL :: INMEMO
      INTEGER(KIND=JPIM) :: N
      INTEGER(KIND=JPIM) :: IMPRES
      INTEGER(KIND=JPIM) :: IO
      INTEGER(KIND=JPIM) :: MODE
      INTEGER(KIND=JPIM) :: NITER
      INTEGER(KIND=JPIM) :: NSIM
      INTEGER(KIND=JPIM) :: M
      INTEGER(KIND=JPIM) :: JMIN
      INTEGER(KIND=JPIM) :: JMAX
      INTEGER(KIND=JPIM) :: IPRINT
      REAL(KIND=JPRB) :: F
      REAL(KIND=JPRB) :: DXMIN
      REAL(KIND=JPRB) :: DF1
      REAL(KIND=JPRB) :: EPSG
      REAL(KIND=JPRB) :: ALPHA(M)
      EXTERNAL :: SIMUL,PROSCA,CTONB,CTCAB
      TYPE(GEOMETRY),    INTENT(INOUT) :: YDGEOMETRY
      TYPE(FIELDS),      INTENT(INOUT) :: YDFIELDS
      TYPE(MTRAJ),       INTENT(INOUT) :: YDMTRAJ
      TYPE(CLASS_VARBC), INTENT(INOUT) :: YDVARBC
      TYPE (CONTROL_VECTOR) :: X
      TYPE (CONTROL_VECTOR) :: G
      TYPE (CONTROL_VECTOR) :: D
      TYPE (CONTROL_VECTOR) :: GG
      TYPE (CONTROL_VECTOR) :: DIAG
      TYPE (CONTROL_VECTOR) :: AUX
      TYPE (CONTROL_VECTOR),DIMENSION(M) :: YBAR
      TYPE (CONTROL_VECTOR),DIMENSION(M) :: SBAR
!
!         variables locales
!
      LOGICAL :: SSCALE,COLD,WARM
      INTEGER(KIND=JPIM) :: I,ITMAX,MODERL,ISIM,JCOUR,INDIC
      REAL(KIND=JPRB) :: R1,T,TMIN,TMAX,GNORM,EPS1,FF
      REAL(KIND=JPRB) :: PRECO,PRECOS,YS,DEN,DK,DK1
      REAL(KIND=JPRB) :: PS,PS2,HP0
      TYPE (CONTROL_VECTOR) :: ZTEMP
!
!         parametres
!
      REAL(KIND=JPRB),PARAMETER :: RM1=0.0001_JPRB
      REAL(KIND=JPRB),PARAMETER :: RM2=0.9_JPRB
      REAL(KIND=JPRB),PARAMETER :: PI=3.1415927_JPRB
      REAL(KIND=JPRB) :: RMIN
!
!---- initialisation
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('M1QN3A',0,ZHOOK_HANDLE)
!
      CALL ALLOCATE_CTLVEC(ZTEMP,X)
      RMIN=1.E-20_JPRB
!
      SSCALE=.TRUE.
      IF (MOD(MODE,2).EQ.0) SSCALE=.FALSE.
!
      WARM=.FALSE.
      IF (MODE/2.EQ.1) WARM=.TRUE.
      COLD=.NOT.WARM
!
      ITMAX=NITER
      NITER=0
      ISIM=1
      EPS1=1.0_JPRB
!
      CALL PROSCA (N,G,G,PS)

      GNORM=SQRT(PS)
      IF (IMPRES.GE.1) WRITE (IO,900) F,GNORM
  900 FORMAT (5X,"f         = ",E15.8                                   &
     &       /5X,"norm of g = ",E15.8)
      IF (GNORM.LT.RMIN) THEN
          MODE=2
          IF (IMPRES.GE.1) WRITE (IO,901)
          GOTO 1000
      ENDIF
  901 FORMAT (/" >>> m1qn3a: initial gradient is too small")
!
!     --- initialisation pour dd
!
      IF (COLD) THEN
          JMIN=1
          JMAX=0
      ENDIF
      JCOUR=1
      IF (INMEMO) JCOUR=JMAX
!
!     --- mise a l'echelle de la premiere direction de descente
!
      IF (COLD) THEN
!
!         --- use Fletcher's scaling and initialize diag to 1.
!
          PRECOS=2.0_JPRB*DF1/GNORM**2
          D%DATA  = -G%DATA * PRECOS
          DIAG = 1.0_JPRB
!
          IF (IMPRES.GE.5) WRITE(IO,902) PRECOS
  902     FORMAT (/" m1qn3a: descent direction -g: precon = ",E10.3)
      ELSE
!
!         --- use the matrix stored in [diag and] the (y,s) pairs
!
          IF (SSCALE) THEN
              PS = DOT_PRODUCT (YBAR(JCOUR),YBAR(JCOUR))
              PRECOS=1.0_JPRB/PS
          ENDIF
          D%DATA = -G%DATA
          IF (INMEMO) THEN
              CALL DD (CTONB,CTCAB,N,SSCALE,M,D,AUX,JMIN,JMAX,          &
     &                 PRECOS,DIAG,ALPHA,YBAR,SBAR)
          ELSE
              CALL DDS (CTONB,CTCAB,N,SSCALE,M,D,AUX,JMIN,JMAX,         &
     &                  PRECOS,DIAG,ALPHA,YBAR(1),SBAR(1))
          ENDIF
      ENDIF
!
      IF (IMPRES.EQ.3) THEN
          WRITE(IO,903)
          WRITE(IO,904)
      ENDIF
      IF (IMPRES.EQ.4) WRITE(IO,903)
  903 FORMAT (/1X,79("-"))
  904 FORMAT (1X)
!
!     --- initialisation pour mlis0
!
      TMAX=1.E+20_JPRB
      CALL PROSCA (N,D,G,HP0)
      IF (HP0.GE.0.0) THEN
          MODE=7
          IF (IMPRES.GE.1) WRITE (IO,905) NITER,HP0
          GOTO 1000
      ENDIF
  905 FORMAT (/" >>> m1qn3 (iteration ",I2,"): "                        &
     &        /5X," the search direction d is not a ",                  &
     &         "descent direction: (g,d) = ",D12.5)
!
!     --- compute the angle (-g,d)
!
      IF (WARM.AND.IMPRES.GE.5) THEN
          CALL PROSCA (N,G,G,PS)
          PS=SQRT(PS)
          CALL PROSCA (N,D,D,PS2)
          PS2=SQRT(PS2)
          PS=HP0/PS/PS2
          PS=MIN(-PS,1.0_JPRB)
          PS=ACOS(PS)
          R1=PS
          R1=R1*180.0_JPRB/PI
          WRITE (IO,906) R1
      ENDIF
  906 FORMAT (/" m1qn3: descent direction d: ",                         &
     &        "angle(-g,d) = ",F5.1," degrees")

!---- Initialize run-time predictor

      CALL PREDICT_RUNTIME (0,ITMAX,EPS1,EPSG)
!
!---- Debut de l'iteration. on cherche x(k+1) de la forme x(k) + t*d,
!     avec t > 0. On connait d.
!
!         Debut de la boucle: etiquette 100,
!         Sortie de la boucle: goto 1000.
!

      ITER_LOOP:DO

      NITER=NITER+1
      IF (IMPRES.LT.0) THEN
          IF (MOD(NITER,-IMPRES).EQ.0) THEN
              INDIC=1
              CALL SIMUL (YDGEOMETRY,YDFIELDS,YDMTRAJ,INDIC,            &
     &                    N,X,F,G,NITER,YDVARBC)
              CYCLE
          ENDIF
      ENDIF
      IF (IMPRES.GE.5) WRITE(IO,903)
      IF (IMPRES.GE.4) WRITE(IO,904)
      IF (IMPRES.GE.3) WRITE (IO,910) NITER,ISIM,F,HP0
  910 FORMAT (" m1qn3: iter ",I3,", simul ",I3,                         &
     &        ", f=",E15.8,", h'(0)=",D12.5)

      GG = G
      FF=F
!
!     --- recherche lineaire et nouveau point x(k+1)
!
      IF (IMPRES.GE.5) WRITE (IO,911)
  911 FORMAT (/" m1qn3: line search")
!
!         --- calcul de tmin
!
      ZTEMP%DATA=ABS(D%DATA)
      TMIN=MAXVAL(ZTEMP)
      TMIN=DXMIN/TMIN
      T=1.
      R1=HP0
!
      CALL MLIS0 (YDGEOMETRY,YDFIELDS,YDMTRAJ,N,SIMUL,PROSCA,           &
     &            X,F,R1,T,TMIN,TMAX,D,G,RM2,RM1,                       &
     &            IMPRES,IO,MODERL,ISIM,NSIM,AUX,NITER,ZTEMP,YDVARBC)
!
!         --- mlis0 renvoie les nouvelles valeurs de x, f et g
!
      IF (MODERL.NE.0) THEN
          IF (MODERL.LT.0) THEN
!
!             --- calcul impossible
!                 t, g: ou les calculs sont impossibles
!                 x, f: ceux du t_gauche (donc f <= ff)
!
              MODE=MODERL
          ELSEIF (MODERL.EQ.1) THEN
!
!             --- descente bloquee sur tmax
!                 [sortie rare (!!) d'apres le code de mlis0]
!
              MODE=3
              IF (IMPRES.GE.1) WRITE(IO,912) NITER
  912         FORMAT (/" >>> m1qn3 (iteration ",I3,                     &
     &                "): line search blocked on tmax: ",               &
     &                "decrease the scaling")
          ELSEIF (MODERL.EQ.4) THEN
!
!             --- nsim atteint
!                 x, f: ceux du t_gauche (donc f <= ff)
!
              MODE=5
          ELSEIF (MODERL.EQ.5) THEN
!
!             --- arret demande par l'utilisateur (indic = 0)
!                 x, f: ceux en sortie du simulateur
!
              MODE=0
          ELSEIF (MODERL.EQ.6) THEN
!
!             --- arret sur dxmin ou appel incoherent
!                 x, f: ceux du t_gauche (donc f <= ff)
!
              MODE=6
          ENDIF
          EXIT
      ENDIF
!
! NOTE: stopping tests are now done after having updated the matrix, so
! that update information can be stored in case of a later warm restart
!
!     --- mise a jour de la matrice
!
      IF (M.GT.0) THEN
!
!         --- mise a jour des pointeurs
!
          JMAX=JMAX+1
          IF (JMAX.GT.M) JMAX=JMAX-M
          IF ((COLD.AND.NITER.GT.M).OR.(WARM.AND.JMIN.EQ.JMAX)) THEN
              JMIN=JMIN+1
              IF (JMIN.GT.M) JMIN=JMIN-M
          ENDIF
          IF (INMEMO) JCOUR=JMAX
!
!         --- y, s et (y,s)
!
          ZTEMP%DATA = T * D%DATA
          CALL CTONB (N,ZTEMP,SBAR(JCOUR))

          ZTEMP%DATA = G%DATA - GG%DATA
          CALL CTONB (N,ZTEMP,YBAR(JCOUR))

          IF (IMPRES.GE.5) THEN
              PS = DOT_PRODUCT (SBAR(JCOUR),SBAR(JCOUR))
              DK1=SQRT(PS)
              IF (NITER.GT.1) WRITE (IO,930) DK1/DK
  930         FORMAT (/" m1qn3: convergence rate, s(k)/s(k-1) = ",      &
     &                E12.5)
              DK=DK1
          ENDIF
          PS = DOT_PRODUCT (YBAR(JCOUR),SBAR(JCOUR))
          YS=PS
          IF (YS.LE.0.0_JPRB) THEN
              MODE=7
              IF (IMPRES.GE.1) WRITE (IO,931) NITER,YS
  931         FORMAT (/" >>> m1qn3 (iteration ",I2,                     &
     &                "): the scalar product (y,s) = ",E12.5            &
     &                /27X,"is not positive")
              EXIT
          ENDIF
!
!         --- ybar et sbar
!
          R1=SQRT(1.0_JPRB/YS)
          SBAR(JCOUR)%DATA = R1 * SBAR(JCOUR)%DATA
          YBAR(JCOUR)%DATA = R1 * YBAR(JCOUR)%DATA
          IF (.NOT.INMEMO) STOP 'ystbl:dead code'
!         if (.not.inmemo) call ystbl (.true.,ybar,sbar,n,jmax)
!
!         --- compute the scalar or diagonal preconditioner
!
          IF (IMPRES.GE.5) WRITE(IO,932)
  932     FORMAT (/" m1qn3: matrix update:")
!
!             --- Here is the Oren-Spedicato factor, for scalar scaling
!
          IF (SSCALE) THEN
              PS = DOT_PRODUCT (YBAR(JCOUR),YBAR(JCOUR))
              PRECOS=1.0_JPRB/PS
!
              IF (IMPRES.GE.5) WRITE (IO,933) PRECOS
  933         FORMAT (5X,"Oren-Spedicato factor = ",E10.3)
!
!             --- Scale the diagonal to Rayleigh's ellipsoid.
!                 Initially (niter.eq.1) and for a cold start, this is
!                 equivalent to an Oren-Spedicato scaling of the
!                 identity matrix.
!
          ELSE
              AUX = YBAR(JCOUR)
              PS=0.0_JPRB
              ZTEMP%DATA = AUX%DATA ** 2
              PS = DOT_PRODUCT(DIAG,ZTEMP)
              R1=1.0_JPRB/PS
              IF (IMPRES.GE.5) THEN
                  WRITE (IO,934) R1
  934             FORMAT(5X,"fitting the ellipsoid: factor = ",E10.3)
              ENDIF
              DIAG%DATA = DIAG%DATA * R1
!
!             --- update the diagonal
!                 (gg is used as an auxiliary vector)
!
              GG = SBAR(JCOUR)
              PS=0.0_JPRB
              ZTEMP%DATA = GG%DATA / DIAG%DATA
              PS = DOT_PRODUCT(GG,ZTEMP)
              DEN=PS
              IPRINT=0

!OCL NOALIAS
              DO I=1,DIAG%NSIZEL
                  DIAG%DATA(I)=1./                                      &
     &            (1./DIAG%DATA(I)+AUX%DATA(I)**2                       &
     &           -(GG%DATA(I)/DIAG%DATA(I))**2/DEN)
                  IF (DIAG%DATA(I).LE.0.0_JPRB) THEN
                      DIAG%DATA(I)=RMIN
                      IPRINT=I
                  ENDIF
              ENDDO
              IF(IPRINT.NE.0) THEN
                IF (IMPRES.GE.5) WRITE (IO,935) IPRINT,                 &
     &           DIAG%DATA(IPRINT),RMIN
              ENDIF

  935         FORMAT (/" >>> m1qn3-WARNING: diagonal element ",I8,      &
     &                 " is negative (",E10.3,"), reset to ",E10.3)
!
              IF (IMPRES.GE.5) THEN
                  PS = SUM(DIAG)
                  PS=PS/N
                  PRECO=PS
!
                  ZTEMP%DATA=(DIAG%DATA-PS)**2
                  PS2=SUM(ZTEMP)
                  PS2=SQRT(PS2/N)
                  WRITE (IO,936) PRECO,PS2
  936             FORMAT (5X,"updated diagonal: average value = ",E10.3,&
     &                   ", sqrt(variance) = ",E10.3)
              ENDIF
          ENDIF
      ENDIF
!
!     --- tests d'arret
!
      CALL PROSCA(N,G,G,PS)
      EPS1=PS
      EPS1=SQRT(EPS1)/GNORM
!
      IF (IMPRES.GE.5) WRITE (IO,940) EPS1
  940 FORMAT (/" m1qn3: stopping criterion on g: ",E12.5)
      IF (EPS1.LT.EPSG) THEN
          MODE=1
          EXIT
      ENDIF
      IF (NITER.EQ.ITMAX) THEN
          MODE=4
          IF (IMPRES.GE.1) WRITE (IO,941) NITER
  941     FORMAT (/" >>> m1qn3 (iteration ",I3,                         &
     &            "): maximal number of iterations")
          EXIT
      ENDIF
      IF (ISIM.GT.NSIM) THEN
          MODE=5
          IF (IMPRES.GE.1) WRITE (IO,942) NITER,ISIM
  942     FORMAT (/" >>> m1qn3 (iteration ",I3,"): ",I6,                &
     &            " simulations (maximal number reached)")
          EXIT
      ENDIF
!
!     --- calcul de la nouvelle direction de descente d = - H.g
!
      IF (M.EQ.0) THEN
          PRECO=2.0_JPRB*(FF-F)/(EPS1*GNORM)**2
          D%DATA = -G%DATA * PRECO
      ELSE
          D%DATA = -G%DATA
          IF (INMEMO) THEN
              CALL DD (CTONB,CTCAB,N,SSCALE,M,D,AUX,JMIN,JMAX,          &
     &                 PRECOS,DIAG,ALPHA,YBAR,SBAR)
          ELSE
              CALL DDS (CTONB,CTCAB,N,SSCALE,M,D,AUX,JMIN,JMAX,         &
     &                  PRECOS,DIAG,ALPHA,YBAR(1),SBAR(1))
          ENDIF
      ENDIF
!
!         --- test: la direction d est-elle de descente ?
!             hp0 sera utilise par mlis0
!
      CALL PROSCA (N,D,G,HP0)
      IF (HP0.GE.0.0_JPRB) THEN
          MODE=7
          IF (IMPRES.GE.1) WRITE (IO,905) NITER,HP0
          EXIT
      ENDIF
      IF (IMPRES.GE.5) THEN
          CALL PROSCA (N,G,G,PS)
          PS=SQRT(PS)
          CALL PROSCA (N,D,D,PS2)
          PS2=SQRT(PS2)
          PS=HP0/PS/PS2
          PS=MIN(-PS,1.0_JPRB)
          PS=ACOS(PS)
          R1=PS
          R1=R1*180.0_JPRB/PI
          WRITE (IO,906) R1
      ENDIF

!---- Predict remaining runtime

      CALL PREDICT_RUNTIME (NITER,ITMAX,EPS1,EPSG)
!
!
!---- on poursuit les iterations
!
      ENDDO ITER_LOOP
!
!---- retour
!
 1000 CONTINUE
      NSIM=ISIM
      EPSG=EPS1
      CALL DEALLOCATE_CTLVEC(ZTEMP)
      IF (LHOOK) CALL DR_HOOK('M1QN3A',1,ZHOOK_HANDLE)
      ENDSUBROUTINE M1QN3A
!
!-----------------------------------------------------------------------
!
