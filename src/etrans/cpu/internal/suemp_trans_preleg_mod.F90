! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! 


MODULE SUEMP_TRANS_PRELEG_MOD
CONTAINS
SUBROUTINE SUEMP_TRANS_PRELEG

! Set up distributed environment for the transform package (part 1)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_GEN         ,ONLY : NOUT, NPRINTLEV
USE TPM_DIM         ,ONLY : R
USE TPM_DISTR       ,ONLY : D, NPRTRW, NPRTRV, MYSETW

USE TPMALD_DISTR    ,ONLY : DALD
USE TPMALD_DIM      ,ONLY : RALD
USE TPMALD_FIELDS   ,ONLY : FALD
USE TPMALD_GEO      ,ONLY : GALD

!USE SUWAVEDI_MOD
!USE ABORT_TRANS_MOD

IMPLICIT NONE

        INTEGER(KIND=JPIM) :: JA,JM,JMLOC,JW,JV,ILATPP,IRESTL,IMLOC,IDT,INM,JN,IM,ILAST

        LOGICAL :: LLP1,LLP2

        INTEGER(KIND=JPIM) :: ISPEC(NPRTRW),IMYMS(RALD%NMSMAX+1),IKNTMP(0:RALD%NMSMAX)
        INTEGER(KIND=JPIM) :: IKMTMP(0:R%NSMAX),ISPEC2P
        INTEGER(KIND=JPIM) :: IC(NPRTRW)
        INTEGER(KIND=JPIM) :: IMDIM,IL,IND,IK,IPOS,IKM
        REAL(KIND=JPRB) :: ZLEPDIM
        REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

        !     ------------------------------------------------------------------

        IF (LHOOK) CALL DR_HOOK('SUEMP_TRANS_PRELEG_MOD:SUEMP_TRANS_PRELEG',0,ZHOOK_HANDLE)
        
        IF(.NOT.D%LGRIDONLY) THEN
                
        LLP1 = NPRINTLEV>0
        LLP2 = NPRINTLEV>1
        IF(LLP1) WRITE(NOUT,*) '=== ENTER ROUTINE SUEMP_TRANS_PRELEG ==='

        !*       1.    Initialize partitioning of wave numbers to PEs ! 
        !             ----------------------------------------------

        ALLOCATE(D%NASM0(0:R%NSMAX))
        IF(LLP2)WRITE(NOUT,9) 'D%NASM0 ',SIZE(D%NASM0   ),SHAPE(D%NASM0   )

        ALLOCATE(DALD%NESM0(0:RALD%NMSMAX))
        IF(LLP2)WRITE(NOUT,9) 'DALD%NESM0 ',SIZE(DALD%NESM0   ),SHAPE(DALD%NESM0   )

        ALLOCATE(D%NATM0(0:R%NTMAX))
        IF(LLP2)WRITE(NOUT,9) 'D%NATM0 ',SIZE(D%NATM0   ),SHAPE(D%NATM0   )
        ALLOCATE(D%NUMPP(NPRTRW))
        IF(LLP2)WRITE(NOUT,9) 'D%NUMPP ',SIZE(D%NUMPP   ),SHAPE(D%NUMPP   )
        ALLOCATE(D%NPOSSP(NPRTRW+1))
        IF(LLP2)WRITE(NOUT,9) 'D%NPOSSP',SIZE(D%NPOSSP  ),SHAPE(D%NPOSSP  )

        ALLOCATE(D%NPROCM(0:RALD%NMSMAX))
        IF(LLP2)WRITE(NOUT,9) 'D%NPROCM',SIZE(D%NPROCM  ),SHAPE(D%NPROCM  )

        ALLOCATE(DALD%NPME(0:RALD%NMSMAX))
        IF(LLP2)WRITE(NOUT,9) 'DALD%NPME',SIZE(DALD%NPME),SHAPE(DALD%NPME)
        ALLOCATE(DALD%NCPL2M(0:RALD%NMSMAX))
        IF(LLP2)WRITE(NOUT,9) 'DALD%NCPL2M',SIZE(DALD%NCPL2M),SHAPE(DALD%NCPL2M)
        CALL ELLIPS(R%NSMAX,RALD%NMSMAX,IKNTMP,IKMTMP)
        DALD%NPME(0)=1
        DO JM=1,RALD%NMSMAX
          DALD%NPME(JM)=DALD%NPME(JM-1)+IKNTMP(JM-1)+1
        ENDDO
        DO JM=0,RALD%NMSMAX
          DALD%NCPL2M(JM) = 2*(IKNTMP(JM)+1)
        ENDDO
        ALLOCATE(FALD%RLEPINM(R%NSPEC_G/2))
        IF(LLP2)WRITE(NOUT,9) 'FALD%RLEPINM',SIZE(FALD%RLEPINM),SHAPE(FALD%RLEPINM)
        DO JM=0,RALD%NMSMAX
          DO JN=1,IKNTMP(JM)
            ZLEPDIM=-((REAL(JM,JPRB)**2)*(GALD%EXWN**2)+&
             & (REAL(JN,JPRB)**2)*(GALD%EYWN**2))  
            FALD%RLEPINM(DALD%NPME(JM)+JN)=1./ZLEPDIM
          ENDDO
        ENDDO
        DO JM=1,RALD%NMSMAX
          ZLEPDIM=-(REAL(JM,JPRB)**2)*(GALD%EXWN**2)
          FALD%RLEPINM(DALD%NPME(JM))=1./ZLEPDIM
        ENDDO
        FALD%RLEPINM(DALD%NPME(0))=0.

        D%NUMPP(:) = 0
        ISPEC(:) = 0
        DALD%NESM0(:)=-99

        IMDIM = 0
        IL = 1
        IND = 1
        IK  = 0
        IPOS = 1
        DO JM=0,RALD%NMSMAX
          IK = IK + IND
          IF (IK > NPRTRW) THEN
            IK = NPRTRW
            IND = -1
          ELSEIF (IK < 1) THEN
            IK = 1
            IND = 1
          ENDIF

          IKM =DALD%NCPL2M(JM)/2 -1
          D%NPROCM(JM) = IK
          ISPEC(IK) = ISPEC(IK)+IKM+1
          D%NUMPP(IK) = D%NUMPP(IK)+1
          IF (IK == MYSETW) THEN
            IMDIM = IMDIM + IKM+1
            IMYMS(IL) = JM
            DALD%NESM0(JM) = IPOS
            IPOS = IPOS+(IKM+1)*4
            IL = IL+1
          ENDIF
        ENDDO
        D%NPOSSP(1) = 1
        ISPEC2P = 4*ISPEC(1)
        D%NSPEC2MX = ISPEC2P
        DO JA=2,NPRTRW
          D%NPOSSP(JA) = D%NPOSSP(JA-1)+ISPEC2P
          ISPEC2P = 4*ISPEC(JA)
          D%NSPEC2MX=MAX(D%NSPEC2MX,ISPEC2P)
        ENDDO
        D%NPOSSP(NPRTRW+1) = D%NPOSSP(NPRTRW)+ISPEC2P

        D%NSPEC2 = 4*IMDIM
        D%NSPEC=D%NSPEC2

        D%NUMP  = D%NUMPP (MYSETW)
        ALLOCATE(D%MYMS(D%NUMP))
        IF(LLP2)WRITE(NOUT,9) 'D%MYMS    ',SIZE(D%MYMS   ),SHAPE(D%MYMS   )
        D%MYMS(:) = IMYMS(1:D%NUMP)
        D%NUMTP = D%NUMP 

        ! pointer to the first wave number of a given wave-set in NALLMS array
        ALLOCATE(D%NPTRMS(NPRTRW))
        IF(LLP2)WRITE(NOUT,9) 'D%NPTRMS   ',SIZE(D%NPTRMS   ),SHAPE(D%NPTRMS   )
        D%NPTRMS(:) = 1
        DO JA=2,NPRTRW
          D%NPTRMS(JA) = D%NPTRMS(JA-1)+D%NUMPP(JA-1)
        ENDDO
        !  D%NALLMS :  wave numbers for all wave-set concatenated together to give all
        !            wave numbers in wave-set order.
        ALLOCATE(D%NALLMS(RALD%NMSMAX+1))
        IF(LLP2)WRITE(NOUT,9) 'D%NALLMS   ',SIZE(D%NALLMS   ),SHAPE(D%NALLMS   )
        IC(:) = 0
        DO JM=0,RALD%NMSMAX
          D%NALLMS(IC(D%NPROCM(JM))+D%NPTRMS(D%NPROCM(JM))) = JM
          IC(D%NPROCM(JM)) = IC(D%NPROCM(JM))+1
        ENDDO
        ALLOCATE(D%NDIM0G(0:RALD%NMSMAX))
        IF(LLP2)WRITE(NOUT,9) 'D%NDIM0G   ',SIZE(D%NDIM0G   ),SHAPE(D%NDIM0G   )
        IPOS = 1
        DO JA=1,NPRTRW
          DO JMLOC=1,D%NUMPP(JA)
            IM = D%NALLMS(D%NPTRMS(JA)+JMLOC-1)
            D%NDIM0G(IM) = IPOS
            IPOS = IPOS+2*DALD%NCPL2M(IM)
          ENDDO
        ENDDO

ALLOCATE(D%NLATLS(NPRTRW,NPRTRV))
IF(LLP2)WRITE(NOUT,9) 'D%NLATLS',SIZE(D%NLATLS   ),SHAPE(D%NLATLS )
ALLOCATE(D%NLATLE(NPRTRW,NPRTRV))
IF(LLP2)WRITE(NOUT,9) 'D%NLATLE',SIZE(D%NLATLE   ),SHAPE(D%NLATLE )

D%NLATLS(:,:) = 9999
D%NLATLE(:,:) = -1

ILATPP = R%NDGL/NPRTRW
IRESTL  = R%NDGL-NPRTRW*ILATPP
DO JW=1,NPRTRW
  IF (JW > IRESTL) THEN
    D%NLATLS(JW,1) = IRESTL*(ILATPP+1)+(JA-IRESTL-1)*ILATPP+1
    D%NLATLE(JW,1) = D%NLATLS(JW,1)+ILATPP-1
  ELSE
    D%NLATLS(JW,1) = (JA-1)*(ILATPP+1)+1
    D%NLATLE(JW,1) = D%NLATLS(JW,1)+ILATPP
  ENDIF
ENDDO
ILAST=0
DO JW=1,NPRTRW
  ILATPP = (D%NLATLE(JW,1)-D%NLATLS(JW,1)+1)/NPRTRV
  IRESTL  = (D%NLATLE(JW,1)-D%NLATLS(JW,1)+1)-NPRTRV*ILATPP
  DO JV=1,NPRTRV
    IF (JV > IRESTL) THEN
      D%NLATLS(JW,JV) = IRESTL*(ILATPP+1)+(JV-IRESTL-1)*ILATPP+1+ILAST
      D%NLATLE(JW,JV) = D%NLATLS(JW,JV)+ILATPP-1
    ELSE
      D%NLATLS(JW,JV) = (JV-1)*(ILATPP+1)+1+ILAST
      D%NLATLE(JW,JV) = D%NLATLS(JW,JV)+ILATPP
    ENDIF
  ENDDO
  ILAST=D%NLATLE(JW,NPRTRV)
ENDDO
IF (LLP1) THEN
  DO JW=1,NPRTRW
    DO JV=1,NPRTRV
      WRITE(NOUT,'(" JW=",I6," JV=",I6," D%NLATLS=",I6," D%NLATLE=",I6)')&
         & JW,JV,D%NLATLS(JW,JV),D%NLATLE(JW,JV)
    ENDDO
  ENDDO
ENDIF

ALLOCATE(D%NPMT(0:R%NSMAX))
IF(LLP2)WRITE(NOUT,9) 'D%NPMT   ',SIZE(D%NPMT   ),SHAPE(D%NPMT   )
ALLOCATE(D%NPMS(0:R%NSMAX))
IF(LLP2)WRITE(NOUT,9) 'D%NPMS   ',SIZE(D%NPMS   ),SHAPE(D%NPMS   )
ALLOCATE(D%NPMG(0:R%NSMAX))
IF(LLP2)WRITE(NOUT,9) 'D%NPMG   ',SIZE(D%NPMG   ),SHAPE(D%NPMG   )
IDT = R%NTMAX-R%NSMAX
INM = 0
DO JMLOC=1,D%NUMP
  IMLOC = D%MYMS(JMLOC)

  INM = INM+R%NTMAX+2-IMLOC
ENDDO
INM = 0
DO JM=0,R%NSMAX

  INM = INM+R%NTMAX+2-JM
ENDDO

D%NLEI3D = (R%NLEI3-1)/NPRTRW+1

ENDIF

IF (LHOOK) CALL DR_HOOK('SUEMP_TRANS_PRELEG_MOD:SUEMP_TRANS_PRELEG',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------
9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

END SUBROUTINE SUEMP_TRANS_PRELEG
END MODULE SUEMP_TRANS_PRELEG_MOD
