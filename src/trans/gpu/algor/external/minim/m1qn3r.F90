      SUBROUTINE M1QN3R (YDGEOMETRY,                                    &
     &                   SIMULR,PROSCAR,CTONBR,CTCABR,N,X,F,G,DXMIN,DF1,&
     &                   EPSG,IMPRES,IO,MODE,NITER,NSIM,IZ,RZ,NRZ)
!AUTOPROMOTE

      USE GEOMETRY_MOD , ONLY : GEOMETRY
      USE PARKIND1, ONLY : JPIM, JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!----
!
!     M1QN3, Version 2.0b, December 1993
!     Jean Charles Gilbert, Claude Lemarechal, INRIA.
!
!     M1qn3 has two running modes: the SID (Scalar Initial Scaling) mode
!     and the DIS (Diagonal Initial Scaling) mode. Both do not require
!     the same amount of storage, the same subroutines, ...
!     In the description below, items that differ in the DIS mode with
!     respect to the SIS mode are given in brakets.
!
!     Use the following subroutines:
!         M1QN3A
!         DD, DDS
!         MLIS0 + ECUBE (Dec 88)
!         MUPDTS, YSTBL.
!
!     The following routines are proposed to the user in case the
!     Euclidean scalar product is used:
!         EUCLID, CTONBE, CTCABE.
!
!     La sous-routine M1QN3 est une interface entre le programme
!     appelant et la sous-routine M1QN3A, le minimiseur proprement dit.
!
!     Le module PROSCA est sense realiser le produit scalaire de deux
!     vecteurs de Rn; le module CTONB est sense realiser le changement
!     de coordonnees correspondant au changement de bases: base
!     euclidienne -> base orthonormale (pour le produit scalaire
!     PROSCA); le module CTBAB fait la transformation inverse: base
!     orthonormale -> base euclidienne.
!
!     Iz is an integer working zone for M1QN3A, its dimension is 5.
!     It is formed of 5 scalars that are set by the optimizer:
!         - the dimension of the problem,
!         - a identifier of the scaling mode,
!         - the number of updates,
!         - two pointers.
!
!     Rz est la zone de travail pour M1QN3A, de dimension nrz.
!     Elle est subdivisee en
!         3 [ou 4] vecteurs de dimension n: d,gg,[diag,]aux
!         m scalaires: alpha
!         m vecteurs de dimension n: ybar
!         m vecteurs de dimension n: sbar
!
!     m est alors le plus grand entier tel que
!         m*(2*n+1)+3*n .le. nrz [m*(2*n+1)+4*n .le. nrz)]
!     soit m := (nrz-3*n) / (2*n+1) [m := (nrz-4*n) / (2*n+1)].
!     Il faut avoir m >= 1, donc nrz >= 5n+1 [nrz >= 6n+1].
!
!     A chaque iteration la metrique est formee a partir d'un multiple
!     de l'identite [d'une matrice diagonale] D qui est mise a jour m
!     fois par la formule de BFGS en utilisant les m couples {y,s} les
!     plus recents.
!
!----
!         arguments
!
      IMPLICIT NONE

      TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
      INTEGER(KIND=JPIM) :: N
      INTEGER(KIND=JPIM) :: IMPRES
      INTEGER(KIND=JPIM) :: IO
      INTEGER(KIND=JPIM) :: MODE
      INTEGER(KIND=JPIM) :: NITER
      INTEGER(KIND=JPIM) :: NSIM
      INTEGER(KIND=JPIM) :: IZ(5)
      INTEGER(KIND=JPIM) :: NRZ
      REAL(KIND=JPRB) :: X(N)
      REAL(KIND=JPRB) :: F
      REAL(KIND=JPRB) :: G(N)
      REAL(KIND=JPRB) :: DXMIN
      REAL(KIND=JPRB) :: DF1
      REAL(KIND=JPRB) :: EPSG
      REAL(KIND=JPRB) :: RZ(NRZ)
      EXTERNAL ::  SIMULR,PROSCAR,CTONBR,CTCABR
!
!         variables locales
!
      LOGICAL :: INMEMO,SSCALE
      INTEGER(KIND=JPIM) :: NTRAVU,ID,IGG,IDIAG,IAUX,IALPHA,IYBAR,ISBAR
      INTEGER(KIND=JPIM) :: M,MMEMO
      REAL(KIND=JPRB) :: R1,R2
      REAL(KIND=JPRB) :: PS
!
!---- impressions initiales et controle des arguments
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('M1QN3R',0,ZHOOK_HANDLE)
      ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM)

      IF (IMPRES.GE.1)                                                  &
     &    WRITE (IO,900) N,DXMIN,DF1,EPSG,NITER,NSIM,IMPRES
900   FORMAT (/" M1QN3 (Version 2.0b, December 1993): entry point"/     &
     &    5X,"dimension of the problem (n):",I9/                        &
     &    5X,"absolute precision on x (dxmin):",E9.2/                   &
     &    5X,"expected decrease for f (df1):",E9.2/                     &
     &    5X,"relative precision on g (epsg):",E9.2/                    &
     &    5X,"maximal number of iterations (niter):",I6/                &
     &    5X,"maximal number of simulations (nsim):",I6/                &
     &    5X,"printing level (impres):",I4)
      IF (N.LE.0.OR.NITER.LE.0.OR.NSIM.LE.0.OR.DXMIN.LE.0.0_JPRB        &
     &    .OR.EPSG.LE.0.0_JPRB                                          &
     &    .OR.EPSG.GT.1.0_JPRB.OR.MODE.LT.0.OR.MODE.GT.3) THEN
          MODE=2
          IF (IMPRES.GE.1) WRITE (IO,901)
901       FORMAT (/" >>> m1qn3r: inconsistent call")
          IF (LHOOK) CALL DR_HOOK('M1QN3R',1,ZHOOK_HANDLE)
          RETURN
      ENDIF
!
!---- what method
!
      IF (MOD(MODE,2).EQ.0) THEN
          IF (IMPRES.GE.1) WRITE (IO,920)
  920     FORMAT (/" m1qn3r: Diagonal Initial Scaling mode")
          SSCALE=.FALSE.
      ELSE
          IF (IMPRES.GE.1) WRITE (IO,921)
  921     FORMAT (/" m1qn3r: Scalar Initial Scaling mode")
          SSCALE=.TRUE.
      ENDIF
!
      IF ((NRZ.LT.5*N+1).OR.((.NOT.SSCALE).AND.(NRZ.LT.6*N+1))) THEN
          MODE=2
          IF (IMPRES.GE.1) WRITE (IO,902)
902       FORMAT (/" >>> m1qn3r: not enough memory allocated")
          IF (LHOOK) CALL DR_HOOK('M1QN3R',1,ZHOOK_HANDLE)
          RETURN
      ENDIF
!
!---- Compute m
!
      CALL MUPDTSR (SSCALE,INMEMO,N,M,NRZ)
!
!     --- Check the value of m (if (y,s) pairs in core, m will be >= 1)
!
      IF (M.LT.1) THEN
          MODE=2
          IF (IMPRES.GE.1) WRITE (IO,9020)
 9020     FORMAT (/" >>> m1qn3r: m is set too small in mupdts")
          IF (LHOOK) CALL DR_HOOK('M1QN3R',1,ZHOOK_HANDLE)
          RETURN
      ENDIF
!
!     --- mmemo = number of (y,s) pairs in core memory
!
      MMEMO=1
      IF (INMEMO) MMEMO=M
!
      NTRAVU=2*(2+MMEMO)*N+M
      IF (SSCALE) NTRAVU=NTRAVU-N
      IF (IMPRES.GE.1) WRITE (IO,903) NRZ,NTRAVU,M
903   FORMAT (/5X,"allocated memory (nrz) :",I9/                        &
     &         5X,"used memory :           ",I9/                        &
     &         5X,"number of updates :     ",I9)
      IF (NRZ.LT.NTRAVU) THEN
          MODE=2
          IF (IMPRES.GE.1) WRITE (IO,902)
          IF (LHOOK) CALL DR_HOOK('M1QN3R',1,ZHOOK_HANDLE)
          RETURN
      ENDIF
!
      IF (IMPRES.GE.1) THEN
          IF (INMEMO) THEN
              WRITE (IO,907)
          ELSE
              WRITE (IO,908)
          ENDIF
      ENDIF
907   FORMAT (5X,"(y,s) pairs are stored in core memory")
908   FORMAT (5X,"(y,s) pairs are stored by the user")
!
!---- cold start or warm restart ?
!     check iz: iz(1)=n, iz(2)=(0 if DIS, 1 if SIS),
!               iz(3)=m, iz(4)=jmin, iz(5)=jmax
!
      IF (MODE/2.EQ.0) THEN
          IF (IMPRES.GE.1) WRITE (IO,922)
      ELSE
          IAUX=0
          IF (SSCALE) IAUX=1
          IF (IZ(1).NE.N.OR.IZ(2).NE.IAUX.OR.IZ(3).NE.M.OR.IZ(4).LT.1   &
     &        .OR.IZ(5).LT.1.OR.IZ(4).GT.IZ(3).OR.IZ(5).GT.IZ(3)) THEN
              MODE=2
              IF (IMPRES.GE.1) WRITE (IO,923)
              IF (LHOOK) CALL DR_HOOK('M1QN3R',1,ZHOOK_HANDLE)
              RETURN
          ENDIF
          IF (IMPRES.GE.1) WRITE (IO,924)
      ENDIF
  922 FORMAT (/" m1qn3r: cold start"/)
  923 FORMAT (/" >>> m1qn3r: inconsistent iz for a warm restart")
  924 FORMAT (/" m1qn3r: warm restart"/)
      IZ(1)=N
      IZ(2)=0
      IF (SSCALE) IZ(2)=1
      IZ(3)=M
!
!---- split the working zone rz
!
      IDIAG=1
      IYBAR=IDIAG+N
      IF (SSCALE) IYBAR=1
      ISBAR=IYBAR+N*MMEMO
      ID=ISBAR+N*MMEMO
      IGG=ID+N
      IAUX=IGG+N
      IALPHA=IAUX+N
!
!---- call the optimization code
!
      CALL M1QN3AR (YDGEOMETRY,                                         &
     &              SIMULR,PROSCAR,CTONBR,CTCABR,N,X,F,G,DXMIN,DF1,EPSG,&
     &              IMPRES,IO,MODE,NITER,NSIM,INMEMO,IZ(3),IZ(4),IZ(5), &
     &              RZ(ID),RZ(IGG),RZ(IDIAG),RZ(IAUX),RZ(IALPHA),       &
     &              RZ(IYBAR),RZ(ISBAR))
!
!---- impressions finales
!
      IF (IMPRES.GE.1) WRITE (IO,905) MODE,NITER,NSIM,EPSG
905   FORMAT (/,1X,79("-")/                                             &
     &        /" m1qn3r: output mode is ",I2                            &
     &        /5X,"number of iterations: ",I4                           &
     &        /5X,"number of simulations: ",I6                          &
     &        /5X,"realized relative precision on g: ",E9.2)
      CALL PROSCAR (YDDIM,N,X,X,PS)
      R1=SQRT(PS)
      CALL PROSCAR (YDDIM,N,G,G,PS)
      R2=SQRT(PS)
      IF (IMPRES.GE.1) WRITE (IO,906) R1,F,R2
906   FORMAT (5X,"norm of x = ",E15.8                                   &
     &       /5X,"f         = ",E15.8                                   &
     &       /5X,"norm of g = ",E15.8)

      END ASSOCIATE
      IF (LHOOK) CALL DR_HOOK('M1QN3R',1,ZHOOK_HANDLE)
      ENDSUBROUTINE M1QN3R
!
!-----------------------------------------------------------------------
