      SUBROUTINE M1QN3AR (YDGEOMETRY,                                   &
     &                    SIMULR,PROSCAR,CTONBR,CTCABR,N,X,F,G,DXMIN,   &
     &                    DF1,EPSG,IMPRES,IO,MODE,NITER,NSIM,INMEMO,M,  &
     &                    JMIN,JMAX,D,GG,DIAG,AUX,ALPHA,YBAR,SBAR)
!AUTOPROMOTE

      USE GEOMETRY_MOD , ONLY : GEOMETRY
      USE PARKIND1, ONLY : JPIM, JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!----
!
!     Code d'optimisation proprement dit.
!
!----
!
!         arguments
!
      IMPLICIT NONE

      TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
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
      REAL(KIND=JPRB) :: X(N)
      REAL(KIND=JPRB) :: F
      REAL(KIND=JPRB) :: G(N)
      REAL(KIND=JPRB) :: DXMIN
      REAL(KIND=JPRB) :: DF1
      REAL(KIND=JPRB) :: EPSG
      REAL(KIND=JPRB) :: D(N)
      REAL(KIND=JPRB) :: GG(N)
      REAL(KIND=JPRB) :: DIAG(N)
      REAL(KIND=JPRB) :: AUX(N)
      REAL(KIND=JPRB) :: ALPHA(M)
      REAL(KIND=JPRB) :: YBAR(N,1)
      REAL(KIND=JPRB) :: SBAR(N,1)
      EXTERNAL :: SIMULR,PROSCAR,CTONBR,CTCABR
!
!         variables locales
!
      LOGICAL :: SSCALE,COLD,WARM
      INTEGER(KIND=JPIM) :: I,ITMAX,MODERL,ISIM,JCOUR,INDIC
      REAL(KIND=JPRB) :: R1,T,TMIN,TMAX,GNORM,EPS1,FF
      REAL(KIND=JPRB) :: PRECO,PRECOS,YS,DEN,DK,DK1
      REAL(KIND=JPRB) :: PS,PS2,HP0
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
      IF (LHOOK) CALL DR_HOOK('M1QN3AR',0,ZHOOK_HANDLE)
      ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM)
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
      EPS1=1.
!
      CALL PROSCAR (YDDIM,N,G,G,PS)
      GNORM=SQRT(PS)
      IF (IMPRES.GE.1) WRITE (IO,900) F,GNORM
  900 FORMAT (5X,"f         = ",E15.8                                   &
     &       /5X,"norm of g = ",E15.8)
      IF (GNORM.LT.RMIN) THEN
          MODE=2
          IF (IMPRES.GE.1) WRITE (IO,901)
          GOTO 1000
      ENDIF
  901 FORMAT (/" >>> m1qn3ar: initial gradient is too small")
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
          PRECOS=2.*DF1/GNORM**2
          DO 10 I=1,N
              D(I)=-G(I)*PRECOS
              DIAG(I)=1.
   10     CONTINUE
          IF (IMPRES.GE.5) WRITE(IO,902) PRECOS
  902     FORMAT (/" m1qn3ar: descent direction -g: precon = ",E10.3)
      ELSE
!
!         --- use the matrix stored in [diag and] the (y,s) pairs
!
          IF (SSCALE) THEN
              CALL PROSCAR (YDDIM,N,YBAR(1,JCOUR),YBAR(1,JCOUR),PS)
              PRECOS=1./PS
          ENDIF
          DO 11 I=1,N
              D(I)=-G(I)
  11      CONTINUE
          IF (INMEMO) THEN
              CALL DDR(YDDIM,PROSCAR,CTONBR,CTCABR,N,SSCALE,M,D,AUX,    &
     &                 JMIN,JMAX,PRECOS,DIAG,ALPHA,YBAR,SBAR)
          ELSE
              CALL DDSR(YDDIM,PROSCAR,CTONBR,CTCABR,N,SSCALE,M,D,AUX,   &
     &                  JMIN,JMAX,PRECOS,DIAG,ALPHA,YBAR,SBAR)
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
      CALL PROSCAR (YDDIM,N,D,G,HP0)
      IF (HP0.GE.0.0_JPRB) THEN
          MODE=7
          IF (IMPRES.GE.1) WRITE (IO,905) NITER,HP0
          GOTO 1000
      ENDIF
  905 FORMAT (/" >>> m1qn3r (iteration ",I2,"): "                       &
     &        /5X," the search direction d is not a ",                  &
     &         "descent direction: (g,d) = ",D12.5)
!
!     --- compute the angle (-g,d)
!
      IF (WARM.AND.IMPRES.GE.5) THEN
          CALL PROSCAR (YDDIM,N,G,G,PS)
          PS=SQRT(PS)
          CALL PROSCAR (YDDIM,N,D,D,PS2)
          PS2=SQRT(PS2)
          PS=HP0/PS/PS2
          PS=MIN(-PS,1.0_JPRB)
          PS=ACOS(PS)
          R1=PS*180.0_JPRB/PI
          WRITE (IO,906) R1
      ENDIF
  906 FORMAT (/" m1qn3r: descent direction d: ",                        &
     &        "angle(-g,d) = ",F5.1," degrees")
!
!---- Debut de l'iteration. on cherche x(k+1) de la forme x(k) + t*d,
!     avec t > 0. On connait d.
!
!         Debut de la boucle: etiquette 100,
!         Sortie de la boucle: goto 1000.
!
100   NITER=NITER+1
      IF (IMPRES.LT.0) THEN
          IF (MOD(NITER,-IMPRES).EQ.0) THEN
              INDIC=1
              CALL SIMULR (YDGEOMETRY,INDIC,N,X,F,G,NITER)
              GOTO 100
          ENDIF
      ENDIF
      IF (IMPRES.GE.5) WRITE(IO,903)
      IF (IMPRES.GE.4) WRITE(IO,904)
      IF (IMPRES.GE.3) WRITE (IO,910) NITER,ISIM,F,HP0
  910 FORMAT (" m1qn3r: iter ",I3,", simulr ",I3,                       &
     &        ", f=",E15.8,", h'(0)=",D12.5)
      DO 101 I=1,N
          GG(I)=G(I)
101   CONTINUE
      FF=F
!
!     --- recherche lineaire et nouveau point x(k+1)
!
      IF (IMPRES.GE.5) WRITE (IO,911)
  911 FORMAT (/" m1qn3r: line search")
!
!         --- calcul de tmin
!
      TMIN=0.
      DO 200 I=1,N
          TMIN=MAX(TMIN,ABS(D(I)))
200   CONTINUE
      TMIN=DXMIN/TMIN
      T=1.
      R1=HP0
!
      CALL MLIS0R (YDGEOMETRY,N,SIMULR,PROSCAR,X,F,R1,T,TMIN,TMAX,D,G,   &
     &             RM2,RM1,IMPRES,IO,MODERL,ISIM,NSIM,AUX,NITER)
!
!         --- mlis0r renvoie les nouvelles valeurs de x, f et g
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
  912         FORMAT (/" >>> m1qn3r (iteration ",I3,                    &
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
          GOTO 1000
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
          DO 400 I=1,N
              SBAR(I,JCOUR)=T*D(I)
              YBAR(I,JCOUR)=G(I)-GG(I)
400       CONTINUE
          IF (IMPRES.GE.5) THEN
              CALL PROSCAR (YDDIM,N,SBAR(1,JCOUR),SBAR(1,JCOUR),PS)
              DK1=SQRT(PS)
              IF (NITER.GT.1) WRITE (IO,930) DK1/DK
  930         FORMAT (/" m1qn3r: convergence rate, s(k)/s(k-1) = ",     &
     &                E12.5)
              DK=DK1
          ENDIF
          CALL PROSCAR (YDDIM,N,YBAR(1,JCOUR),SBAR(1,JCOUR),PS)
          YS=PS
          IF (YS.LE.0.) THEN
              MODE=7
              IF (IMPRES.GE.1) WRITE (IO,931) NITER,YS
  931         FORMAT (/" >>> m1qn3r (iteration ",I2,                    &
     &                "): the scalar product (y,s) = ",E12.5            &
     &                /27X,"is not positive")
              GOTO 1000
          ENDIF
!
!         --- ybar et sbar
!
          R1=SQRT(1./YS)
          DO 410 I=1,N
              SBAR(I,JCOUR)=R1*SBAR(I,JCOUR)
              YBAR(I,JCOUR)=R1*YBAR(I,JCOUR)
  410     CONTINUE
          IF (.NOT.INMEMO) CALL YSTBLR (.TRUE.,YBAR,SBAR,N,JMAX)
!
!         --- compute the scalar or diagonal preconditioner
!
          IF (IMPRES.GE.5) WRITE(IO,932)
  932     FORMAT (/" m1qn3r: matrix update:")
!
!             --- Here is the Oren-Spedicato factor, for scalar scaling
!
          IF (SSCALE) THEN
              CALL PROSCAR (YDDIM,N,YBAR(1,JCOUR),YBAR(1,JCOUR),PS)
              PRECOS=1./PS
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
              CALL CTONBR (N,YBAR(1,JCOUR),AUX)
              PS=0.0
              DO 420 I=1,N
                  PS=PS+(DIAG(I)*AUX(I)*AUX(I))
  420         CONTINUE
              R1=1.0/PS
              IF (IMPRES.GE.5) THEN
                  WRITE (IO,934) R1
  934             FORMAT(5X,"fitting the ellipsoid: factor = ",E10.3)
              ENDIF
              DO 421 I=1,N
                  DIAG(I)=DIAG(I)*R1
  421         CONTINUE
!
!             --- update the diagonal
!                 (gg is used as an auxiliary vector)
!
              CALL CTONBR (N,SBAR(1,JCOUR),GG)
              PS=0.0
              DO 430 I=1,N
                  PS=PS+(GG(I)*GG(I)/DIAG(I))
  430         CONTINUE
              DEN=PS
              IPRINT=0
              DO 431 I=1,N
                  DIAG(I)=1.0_JPRB/                                     &
     &               (1.0_JPRB/DIAG(I)+AUX(I)**2-(GG(I)/DIAG(I))**2/DEN)
                  IF (DIAG(I).LE.0.0_JPRB) THEN
!EA                      if (impres.ge.5) write (io,935) i,diag(i),rmin
                      DIAG(I)=RMIN
                      IPRINT=I
                  ENDIF
  431         CONTINUE
              IF(IPRINT.NE.0) THEN
                IF (IMPRES.GE.5) WRITE (IO,935) IPRINT,DIAG(IPRINT),RMIN
              ENDIF
  935         FORMAT (/" >>> m1qn3r-WARNING: diagonal element ",I8,     &
     &                 " is negative (",E10.3,"), reset to ",E10.3)
!
              IF (IMPRES.GE.5) THEN
                  PS=0.
                  DO 440 I=1,N
                      PS=PS+DIAG(I)
  440             CONTINUE
                  PS=PS/N
                  PRECO=PS
!
                  PS2=0.
                  DO 441 I=1,N
                      PS2=PS2+(DIAG(I)-PS)**2
  441             CONTINUE
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
      CALL PROSCAR(YDDIM,N,G,G,PS)
      EPS1=SQRT(PS)/GNORM
!
      IF (IMPRES.GE.5) WRITE (IO,940) EPS1
  940 FORMAT (/" m1qn3r: stopping criterion on g: ",E12.5)
      IF (EPS1.LT.EPSG) THEN
          MODE=1
          GOTO 1000
      ENDIF
      IF (NITER.EQ.ITMAX) THEN
          MODE=4
          IF (IMPRES.GE.1) WRITE (IO,941) NITER
  941     FORMAT (/" >>> m1qn3r (iteration ",I3,                        &
     &            "): maximal number of iterations")
          GOTO 1000
      ENDIF
      IF (ISIM.GT.NSIM) THEN
          MODE=5
          IF (IMPRES.GE.1) WRITE (IO,942) NITER,ISIM
  942     FORMAT (/" >>> m1qn3r (iteration ",I3,"): ",I6,               &
     &            " simulations (maximal number reached)")
          GOTO 1000
      ENDIF
!
!     --- calcul de la nouvelle direction de descente d = - H.g
!
      IF (M.EQ.0) THEN
          PRECO=2.0_JPRB*(FF-F)/(EPS1*GNORM)**2
          DO 500 I=1,N
              D(I)=-G(I)*PRECO
  500     CONTINUE
      ELSE
          DO 510 I=1,N
              D(I)=-G(I)
  510     CONTINUE
          IF (INMEMO) THEN
              CALL DDR(YDDIM,PROSCAR,CTONBR,CTCABR,N,SSCALE,M,D,AUX,    &
     &                 JMIN,JMAX,PRECOS,DIAG,ALPHA,YBAR,SBAR)
          ELSE
              CALL DDSR(YDDIM,PROSCAR,CTONBR,CTCABR,N,SSCALE,M,D,AUX,   &
     &                  JMIN,JMAX,PRECOS,DIAG,ALPHA,YBAR,SBAR)
          ENDIF
      ENDIF
!
!         --- test: la direction d est-elle de descente ?
!             hp0 sera utilise par mlis0
!
      CALL PROSCAR (YDDIM,N,D,G,HP0)
      IF (HP0.GE.0.0) THEN
          MODE=7
          IF (IMPRES.GE.1) WRITE (IO,905) NITER,HP0
          GOTO 1000
      ENDIF
      IF (IMPRES.GE.5) THEN
          CALL PROSCAR (YDDIM,N,G,G,PS)
          PS=SQRT(PS)
          CALL PROSCAR (YDDIM,N,D,D,PS2)
          PS2=SQRT(PS2)
          PS=HP0/PS/PS2
          PS=MIN(-PS,1.0_JPRB)
          PS=ACOS(PS)
          R1=PS*180.0_JPRB/PI
          WRITE (IO,906) R1
      ENDIF
!
!---- on poursuit les iterations
!
      GOTO 100
!
!---- retour
!
 1000 CONTINUE
      NSIM=ISIM
      EPSG=EPS1

      END ASSOCIATE
      IF (LHOOK) CALL DR_HOOK('M1QN3AR',1,ZHOOK_HANDLE)
      ENDSUBROUTINE M1QN3AR
!
!-----------------------------------------------------------------------
