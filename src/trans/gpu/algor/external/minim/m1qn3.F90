      SUBROUTINE M1QN3 (YDGEOMETRY,YDFIELDS,YDMTRAJ,SIMUL,PROSCA,CTONB, &
     &                  CTCAB,N,X,F,G,DXMIN,DF1,                        &
     &                  EPSG,IMPRES,IO,MODE,NITER,NSIM,IZ,M,MMEMO,      &
     &                  RZDIAG,YBAR,SBAR,YDVARBC)
!AUTOPROMOTE
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
!    Modifications:
!    --------------
!
!       2001-07-27: M. Fisher - {y,s} pairs are now held with respect to
!                               the Euclidean inner product. This greatly
!                               reduces the number of calls to prosca.
!                               NB: if you manipulate {y,s} outside of
!                               m1qn3, consider uncommenting the (untested)
!                               code before and after the call to m1qn3a.
!
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
      INTEGER(KIND=JPIM) :: N
      INTEGER(KIND=JPIM) :: IMPRES
      INTEGER(KIND=JPIM) :: IO
      INTEGER(KIND=JPIM) :: MODE
      INTEGER(KIND=JPIM) :: NITER
      INTEGER(KIND=JPIM) :: NSIM
      INTEGER(KIND=JPIM) :: IZ(5)
      INTEGER(KIND=JPIM) :: M
      INTEGER(KIND=JPIM) :: MMEMO
      REAL(KIND=JPRB) :: F
      REAL(KIND=JPRB) :: DXMIN
      REAL(KIND=JPRB) :: DF1
      REAL(KIND=JPRB) :: EPSG
      EXTERNAL :: SIMUL,PROSCA,CTONB,CTCAB
      TYPE(GEOMETRY),   INTENT(INOUT) :: YDGEOMETRY
      TYPE(FIELDS),     INTENT(INOUT) :: YDFIELDS
      TYPE(MTRAJ),      INTENT(INOUT) :: YDMTRAJ
      TYPE(CLASS_VARBC),INTENT(INOUT) :: YDVARBC
      TYPE (CONTROL_VECTOR) :: X
      TYPE (CONTROL_VECTOR) :: G
      TYPE (CONTROL_VECTOR),DIMENSION(MMEMO) :: YBAR
      TYPE (CONTROL_VECTOR),DIMENSION(MMEMO) :: SBAR
      TYPE (CONTROL_VECTOR) :: RZDIAG
!
!         variables locales
!
      TYPE (CONTROL_VECTOR) :: RZD,RZGG,RZAUX
      REAL(KIND=JPRB) :: ALPHA(M)
      LOGICAL :: INMEMO
      LOGICAL :: SSCALE
      INTEGER(KIND=JPIM) :: ID,IGG,IDIAG,IAUX,IALPHA,IYBAR,ISBAR,J,JFIN
      REAL(KIND=JPRB) :: R1,R2
      REAL(KIND=JPRB) :: PS
!
!---- impressions initiales et controle des arguments
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('M1QN3',0,ZHOOK_HANDLE)
      IF (IMPRES.GE.1)                                                  &
     &    WRITE (IO,900) N,M,DXMIN,DF1,EPSG,NITER,NSIM,IMPRES
900   FORMAT (/" M1QN3 (Version 2.0b, December 1993): entry point"/     &
     &    5X,"dimension of the problem (n):",I9/                        &
     &    5X,"number of updates (m):",I9/                               &
     &    5X,"absolute precision on x (dxmin):",E9.2/                   &
     &    5X,"expected decrease for f (df1):",E9.2/                     &
     &    5X,"relative precision on g (epsg):",E9.2/                    &
     &    5X,"maximal number of iterations (niter):",I6/                &
     &    5X,"maximal number of simulations (nsim):",I6/                &
     &    5X,"printing level (impres):",I4)
      IF (N.LE.0.OR.NITER.LE.0.OR.NSIM.LE.0.OR.DXMIN.LE.0.0_JPRB        &
     &    .OR.EPSG.LE.0.0_JPRB                                          &
     &    .OR.EPSG.GT.1.0_JPRB .OR.MODE.LT.0.OR.MODE.GT.3) THEN
          MODE=2
          IF (IMPRES.GE.1) WRITE (IO,901)
901       FORMAT (/" >>> m1qn3: inconsistent call")
          IF (LHOOK) CALL DR_HOOK('M1QN3',1,ZHOOK_HANDLE)
          RETURN
      ENDIF
!
!---- what method
!
      IF (MOD(MODE,2).EQ.0) THEN
          IF (IMPRES.GE.1) WRITE (IO,920)
  920     FORMAT (/" m1qn3: Diagonal Initial Scaling mode")
          SSCALE=.FALSE.
      ELSE
          IF (IMPRES.GE.1) WRITE (IO,921)
  921     FORMAT (/" m1qn3: Scalar Initial Scaling mode")
          SSCALE=.TRUE.
      ENDIF
!
!     --- Check the value of m (if (y,s) pairs in core, m will be >= 1)
!
      IF (M.LT.1) THEN
          MODE=2
          IF (IMPRES.GE.1) WRITE (IO,9020)
 9020     FORMAT (/" >>> m1qn3: m is set too small in mupdts")
          IF (LHOOK) CALL DR_HOOK('M1QN3',1,ZHOOK_HANDLE)
          RETURN
      ENDIF
!
!     --- mmemo = number of (y,s) pairs in core memory
!
      CALL ALLOCATE_CTLVEC(RZD,  X)
      CALL ALLOCATE_CTLVEC(RZGG, X)
      CALL ALLOCATE_CTLVEC(RZAUX,X)
!
      IF(M.EQ.MMEMO) THEN
        INMEMO=.TRUE.
      ELSE
        INMEMO=.FALSE.
      ENDIF
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
              IF (LHOOK) CALL DR_HOOK('M1QN3',1,ZHOOK_HANDLE)
              RETURN
          ENDIF
          IF (IMPRES.GE.1) WRITE (IO,924)
      ENDIF
  922 FORMAT (/" m1qn3: cold start"/)
  923 FORMAT (/" >>> m1qn3: inconsistent iz for a warm restart")
  924 FORMAT (/" m1qn3: warm restart"/)
      IZ(1)=N
      IZ(2)=0
      IF (SSCALE) IZ(2)=1
      IZ(3)=M
!
!---- split the working zone rz
!
!      idiag=1
!      iybar=idiag+n
!      if (sscale) iybar=1
!      isbar=iybar+n*mmemo
!      id=isbar+n*mmemo
!      igg=id+n
!      iaux=igg+n
!      ialpha=iaux+n

!=========== uncomment the following if you require m1qn3 to
!=========== input and output the (y,s) pairs w.r.t. the prosca inner
!=========== product. If you don't care how (y,s) are stored externally
!=========== then this transformation is unnecessary expense.
!
!---- transform the (y,s) pairs to Euclidean inner product
!
!     if (mode/2.ne.0) then
!       jfin=jmax
!       if (jfin.lt.jmin) jfin=jmax+iz(3)
!       do j=jfin,jmin,-1
!         jp=j
!         if (jp.gt.iz(3)) jp=jp-iz(3)
!         if (inmemo) then
!           rzaux=sbar(jp)
!           call ctonb (n,rzaux,sbar(jp))
!           rzaux=ybar(jp)
!           call ctonb (n,rzaux,ybar(jp))
!         else
!           call ystbl (.false.,ybar,sbar,n,jp)
!           rzaux=sbar(1)
!           call ctonb (n,rzaux,sbar(1))
!           rzaux=ybar(1)
!           call ctonb (n,rzaux,ybar(1))
!           call ystbl (.true.,ybar,sbar,n,jp)
!         endif
!       enddo
!     endif
!
!---- call the optimization code
!
      CALL M1QN3A (YDGEOMETRY,YDFIELDS,YDMTRAJ,SIMUL,PROSCA,CTONB,      &
     &             CTCAB,N,X,F,G,DXMIN,DF1,EPSG,                        &
     &             IMPRES,IO,MODE,NITER,NSIM,INMEMO,IZ(3),IZ(4),IZ(5),  &
     &             RZD,RZGG,RZDIAG,RZAUX,ALPHA,                         &
     &             YBAR,SBAR,YDVARBC)

!=========== uncomment the following if you require m1qn3 to
!=========== input and output the (y,s) pairs w.r.t. the prosca inner
!=========== product. If you don't care how (y,s) are stored externally
!=========== then this transformation is unnecessary expense.
!
!---- transform the (y,s) pairs back to a prosca inner product
!
!     if (mode/2.ne.0) then
!       jfin=jmax
!       if (jfin.lt.jmin) jfin=jmax+iz(3)
!       do j=jfin,jmin,-1
!         jp=j
!         if (jp.gt.iz(3)) jp=jp-iz(3)
!         if (inmemo) then
!           rzaux=sbar(jp)
!           call ctcab (n,rzaux,sbar(jp))
!           rzaux=ybar(jp)
!           call ctcab (n,rzaux,ybar(jp))
!         else
!           call ystbl (.false.,ybar,sbar,n,jp)
!           rzaux=sbar(1)
!           call ctcab (n,rzaux,sbar(1))
!           rzaux=ybar(1)
!           call ctcab (n,rzaux,ybar(1))
!           call ystbl (.true.,ybar,sbar,n,jp)
!         endif
!       enddo
!     endif
!
!---- impressions finales
!
      IF (IMPRES.GE.1) WRITE (IO,905) MODE,NITER,NSIM,EPSG
905   FORMAT (/,1X,79("-")/                                             &
     &        /" m1qn3: output mode is ",I2                             &
     &        /5X,"number of iterations: ",I4                           &
     &        /5X,"number of simulations: ",I6                          &
     &        /5X,"realized relative precision on g: ",E9.2)
      CALL PROSCA (N,X,X,PS)
      R1=SQRT(PS)
      CALL PROSCA (N,G,G,PS)
      R2=SQRT(PS)
      IF (IMPRES.GE.1) WRITE (IO,906) R1,F,R2
906   FORMAT (5X,"norm of x = ",E15.8                                   &
     &       /5X,"f         = ",E15.8                                   &
     &       /5X,"norm of g = ",E15.8)

      CALL DEALLOCATE_CTLVEC(RZD)
      CALL DEALLOCATE_CTLVEC(RZGG)
      CALL DEALLOCATE_CTLVEC(RZAUX)
      IF (LHOOK) CALL DR_HOOK('M1QN3',1,ZHOOK_HANDLE)
      ENDSUBROUTINE M1QN3
!
!-----------------------------------------------------------------------
!
