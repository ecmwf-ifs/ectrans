MODULE DIR_TRANS_CTLAD_MOD
CONTAINS
SUBROUTINE DIR_TRANS_CTLAD(KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS,KF_UV,KF_SCALARS,&
 & PSPVOR,PSPDIV,PSPSCALAR,KVSETUV,KVSETSC,PGP)

!**** *DIR_TRANS_CTLAD* - Control routine for direct spectral transform-adj.

!     Purpose.
!     --------
!        Control routine for the direct spectral transform

!**   Interface.
!     ----------
!     CALL DIR_TRANS_CTLAD(...)

!     Explicit arguments : 
!     -------------------- 
!     KF_UV_G      - global number of spectral u-v fields
!     KF_SCALARS_G - global number of scalar spectral fields
!     KF_GP        - total number of output gridpoint fields
!     KF_FS        - total number of fields in fourier space
!     KF_UV        - local number of spectral u-v fields
!     KF_SCALARS   - local number of scalar spectral fields
!     PSPVOR(:,:)  - spectral vorticity
!     PSPDIV(:,:)  - spectral divergence
!     PSPSCALAR(:,:) - spectral scalarvalued fields
!     KVSETUV(:)  - indicating which 'b-set' in spectral space owns a 
!                   vor/div field. Equivalant to NBSETLEV in the IFS.
!                   The length of KVSETUV should be the GLOBAL number
!                   of u/v fields which is the dimension of u and v releated
!                   fields in grid-point space. 
!     KVESETSC(:) - indicating which 'b-set' in spectral space owns a
!                   scalar field. As for KVSETUV this argument is required
!                   if the total number of processors is greater than
!                   the number of processors used for distribution in
!                   spectral wave space.  
!     PGP(:,:,:)  - gridpoint fields 

!                  The ordering of the output fields is as follows (all 
!                  parts are optional depending on the input switches):
!
!       u             : KF_UV_G fields
!       v             : KF_UV_G fields
!       scalar fields : KF_SCALARS_G fields 
!   
!     Method.
!     -------

!     Externals.  SHUFFLE     - reshuffle fields for load balancing
!     ----------  FIELD_SPLIT - split fields in NPROMATR packets
!                 LTDIR_CTLAD   - control of Legendre transform
!                 FTDIR_CTLAD   - control of Fourier transform

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 01-01-03

!     ------------------------------------------------------------------

#include "tsmbkind.h"

USE TPM_GEN
USE TPM_TRANS
USE TPM_DISTR

USE SHUFFLE_MOD
USE FIELD_SPLIT_MOD
USE LTDIR_CTLAD_MOD
USE FTDIR_CTLAD_MOD

IMPLICIT NONE

! Declaration of arguments

INTEGER_M, INTENT(IN) :: KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS,KF_UV,KF_SCALARS
REAL_B    ,OPTIONAL, INTENT(INOUT) :: PSPVOR(:,:)
REAL_B    ,OPTIONAL, INTENT(INOUT) :: PSPDIV(:,:)
REAL_B    ,OPTIONAL, INTENT(INOUT) :: PSPSCALAR(:,:)
INTEGER_M ,OPTIONAL, INTENT(IN)  :: KVSETUV(:)
INTEGER_M ,OPTIONAL, INTENT(IN)  :: KVSETSC(:)
REAL_B    ,INTENT(OUT) :: PGP(:,:,:)

! Local variables

INTEGER_M :: IPTRGP(KF_GP),IPTRSPUV(NPROMATR),IPTRSPSC(NPROMATR)
INTEGER_M :: ISHFUV_G(KF_GP),ISHFSC_G(KF_GP)
INTEGER_M :: IVSETUV(KF_GP),IVSETSC(KF_GP)
INTEGER_M :: IBLKS,J,JBLK,ISTF,IENF,IFLDS,ISTUV_G,IENUV_G
INTEGER_M :: IF_UV_G,IF_UV,ISTUV,IF_SCALARS,IF_SCALARS_G,IF_FS,IF_GP
INTEGER_M :: IFLDSUV,JFLD,ISTSC_G,IENSC_G,ISTSC,IENSC,IENUV,IF_GPB


!     ------------------------------------------------------------------

! Perform transform

IF_GPB = 2*KF_UV_G+KF_SCALARS_G
IF(NPROMATR > 0 .AND. IF_GPB > NPROMATR) THEN

  ! Fields to be split into packets

  CALL SHUFFLE(KF_UV_G,KF_SCALARS_G,ISHFUV_G,IVSETUV,ISHFSC_G,IVSETSC,& 
 & KVSETUV,KVSETSC)

  IBLKS=(IF_GPB-1)/NPROMATR+1

  DO JBLK=1,IBLKS
  
    CALL FIELD_SPLIT(JBLK,KF_GP,KF_UV_G,IVSETUV,IVSETSC,&
     & ISTUV_G,IENUV_G,IF_UV_G,ISTSC_G,IENSC_G,IF_SCALARS_G,&
     & ISTUV,IENUV,IF_UV,ISTSC,IENSC,IF_SCALARS)

    IF_FS = 2*IF_UV + IF_SCALARS
    IF_GP = 2*IF_UV_G+IF_SCALARS_G
    DO JFLD=1,IF_UV_G
      IPTRGP(JFLD) = ISHFUV_G(ISTUV_G+JFLD-1)
      IPTRGP(JFLD+IF_UV_G) = KF_UV_G+ISHFUV_G(ISTUV_G+JFLD-1)
    ENDDO
    DO JFLD=1,IF_SCALARS_G
      IPTRGP(JFLD+2*IF_UV_G) = 2*KF_UV_G+ISHFSC_G(ISTSC_G+JFLD-1)
    ENDDO
    DO JFLD=1,IF_UV
      IPTRSPUV(JFLD) = ISTUV+JFLD-1
    ENDDO
    DO JFLD=1,IF_SCALARS
      IPTRSPSC(JFLD) = ISTSC+JFLD-1
    ENDDO

    CALL LTDIR_CTLAD(IF_FS,IF_UV,IF_SCALARS, &
     & PSPVOR=PSPVOR,PSPDIV=PSPDIV,PSPSCALAR=PSPSCALAR,&
     & KFLDPTRUV=IPTRSPUV,KFLDPTRSC=IPTRSPSC)
    IF(IF_UV_G > 0 .AND. IF_SCALARS_G > 0) THEN
      CALL FTDIR_CTLAD(PGP,IF_UV_G,IF_SCALARS_G,IF_GP,IF_FS,&
       & KVSETUV=IVSETUV(ISTUV_G:IENUV_G),&
       & KVSETSC=IVSETSC(ISTSC_G:IENSC_G),KPTRGP=IPTRGP)
    ELSEIF(IF_UV_G > 0) THEN
      CALL FTDIR_CTLAD(PGP,IF_UV_G,IF_SCALARS_G,IF_GP,IF_FS,&
       & KVSETUV=IVSETUV(ISTUV_G:IENUV_G),&
       & KPTRGP=IPTRGP)
    ELSEIF(IF_SCALARS_G > 0) THEN
      CALL FTDIR_CTLAD(PGP,IF_UV_G,IF_SCALARS_G,IF_GP,IF_FS,&
       & KVSETSC=IVSETSC(ISTSC_G:IENSC_G),KPTRGP=IPTRGP)
    ENDIF
  ENDDO

ELSE

  ! No splitting of fields, transform done in one go

  CALL LTDIR_CTLAD(KF_FS,KF_UV,KF_SCALARS, &
   &PSPVOR,PSPDIV,PSPSCALAR)
  CALL FTDIR_CTLAD(PGP,KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS,&
   & KVSETUV,KVSETSC)

ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE DIR_TRANS_CTLAD
END MODULE DIR_TRANS_CTLAD_MOD
