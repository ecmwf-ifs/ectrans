MODULE INV_TRANS_CTLAD_MOD
CONTAINS
SUBROUTINE INV_TRANS_CTLAD(KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS,KF_OUT_LT,&
 & KF_UV,KF_SCALARS,KF_SCDERS,&
 & PSPVOR,PSPDIV,PSPSCALAR,KVSETUV,KVSETSC,PGP)

!**** *INV_TRANS_CTLAD* - Control routine for inverse spectral transform adj.

!     Purpose.
!     --------
!        Control routine for the inverse spectral transform

!**   Interface.
!     ----------
!     CALL INV_TRANS_CTLAD(...)

!     Explicit arguments : 
!     -------------------- 
!     KF_UV_G      - global number of spectral u-v fields
!     KF_SCALARS_G - global number of scalar spectral fields
!     KF_GP        - total number of output gridpoint fields
!     KF_FS        - total number of fields in fourier space
!     KF_OUT_LT    - total number of fields coming out from inverse LT
!     KF_UV        - local number of spectral u-v fields
!     KF_SCALARS   - local number of scalar spectral fields
!     KF_SCDERS    - local number of derivatives of scalar spectral fields
!     PSPVOR(:,:)  - spectral vorticity (input)
!     PSPDIV(:,:)  - spectral divergence (input)
!     PSPSCALAR(:,:) - spectral scalarvalued fields (input)
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
!     PGP(:,:,:)  - gridpoint fields (output)

!                  The ordering of the output fields is as follows (all 
!                  parts are optional depending on the input switches):
!
!       vorticity     : KF_UV_G fields
!       divergence    : KF_UV_G fields
!       u             : KF_UV_G fields
!       v             : KF_UV_G fields
!       scalar fields : KF_SCALARS_G fields 
!       N-S derivative of scalar fields : KF_SCALARS_G fields 
!       E-W derivative of u : KF_UV_G fields 
!       E-W derivative of v : KF_UV_G fields 
!       E-W derivative of scalar fields : KF_SCALARS_G fields 
!   
!     Method.
!     -------

!     Externals.  SHUFFLE     - reshuffle fields for load balancing
!     ----------  FIELD_SPLIT - split fields in NPROMATR packets
!                 LTINV_CTLAD   - control of Legendre transform
!                 FTINV_CTLAD   - control of Fourier transform

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
USE LTINV_CTLAD_MOD
USE FTINV_CTLAD_MOD

IMPLICIT NONE

! Declaration of arguments

INTEGER_M, INTENT(IN) :: KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS,KF_OUT_LT
INTEGER_M, INTENT(IN) :: KF_UV,KF_SCALARS,KF_SCDERS
REAL_B    ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
REAL_B    ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
REAL_B    ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)
INTEGER_M ,OPTIONAL, INTENT(IN) :: KVSETUV(:)
INTEGER_M ,OPTIONAL, INTENT(IN) :: KVSETSC(:)
REAL_B    ,INTENT(IN) :: PGP(:,:,:)

! Local variables

INTEGER_M :: IPTRGP(KF_GP),IPTRSPUV(NPROMATR),IPTRSPSC(NPROMATR)
INTEGER_M :: ISHFUV_G(KF_GP),ISHFSC_G(KF_GP)
INTEGER_M :: IVSETUV(KF_GP),IVSETSC(KF_GP)
INTEGER_M :: IBLKS,J,JBLK,ISTF,IENF,IFLDS,ISTUV_G,IENUV_G
INTEGER_M :: IF_UV_G,IF_UV,ISTUV,IF_SCALARS,IF_SCALARS_G,IF_FS,IF_GP
INTEGER_M :: IFLDSUV,JFLD,ISTSC_G,IENSC_G,ISTSC,IENSC,IENUV,IF_SCDERS,IF_OUT_LT
INTEGER_M :: IOFFD,IOFFU,IOFFV,IOFFUVD,IOFFSC,IOFFSCNS,IOFFSCEW,IOFF,IF_GPB

!     ------------------------------------------------------------------

! Perform transform


IF_GPB = 2*KF_UV_G+KF_SCALARS_G
IF(NPROMATR > 0 .AND. IF_GPB > NPROMATR) THEN

  ! Fields to be split into packets

  CALL SHUFFLE(KF_UV_G,KF_SCALARS_G,ISHFUV_G,IVSETUV,ISHFSC_G,IVSETSC,& 
 & KVSETUV,KVSETSC)

  IBLKS=(IF_GPB-1)/NPROMATR+1

  DO JBLK=1,IBLKS
  
    CALL FIELD_SPLIT(JBLK,IF_GPB,KF_UV_G,IVSETUV,IVSETSC,&
     & ISTUV_G,IENUV_G,IF_UV_G,ISTSC_G,IENSC_G,IF_SCALARS_G,&
     & ISTUV,IENUV,IF_UV,ISTSC,IENSC,IF_SCALARS)

    IF(LSCDERS) THEN
      IF_SCDERS = IF_SCALARS
    ELSE
      IF_SCDERS = 0
    ENDIF
        
    IF_OUT_LT = 2*IF_UV + IF_SCALARS+IF_SCDERS
    IF(LVORGP) THEN
      IF_OUT_LT = IF_OUT_LT+IF_UV
    ENDIF
    IF(LDIVGP) THEN
      IF_OUT_LT = IF_OUT_LT+IF_UV
    ENDIF
    IF_FS = IF_OUT_LT+IF_SCDERS
    IF(LUVDER) THEN
      IF_FS = IF_FS+2*IF_UV
    ENDIF

    IF_GP = 2*IF_UV_G+IF_SCALARS_G
    IOFFD = 0
    IOFFU = 0
    IOFFV = KF_UV_G
    IOFFUVD = 2*KF_UV_G+KF_SCALARS_G
    IOFFSC = 2*KF_UV_G
    IF(LVORGP) THEN
      IF_GP = IF_GP+IF_UV_G
      IOFFD = KF_UV_G
      IOFFU = IOFFU+KF_UV_G
      IOFFV = IOFFV+KF_UV_G
      IOFFUVD =IOFFUVD+KF_UV_G
      IOFFSC = IOFFSC+KF_UV_G
    ENDIF
    IF(LDIVGP) THEN
      IF_GP = IF_GP+IF_UV_G
      IOFFU = IOFFU+KF_UV_G
      IOFFV = IOFFV+KF_UV_G
      IOFFUVD =IOFFUVD+KF_UV_G
      IOFFSC = IOFFSC+KF_UV_G
    ENDIF
    IF(LSCDERS) THEN
      IF_GP  = IF_GP+2*IF_SCALARS_G
      IOFFUVD =IOFFUVD+KF_SCALARS_G
      IOFFSCNS = IOFFSC+KF_SCALARS_G
      IOFFSCEW = IOFFSC+2*KF_SCALARS_G
    ENDIF
    IF(LUVDER) THEN
      IF_GP = IF_GP+2*IF_UV_G
      IOFFSCEW = IOFFSCEW+2*KF_UV_G
    ENDIF

    DO JFLD=1,IF_UV_G
      IOFF = 0
      IF(LVORGP) THEN
        IPTRGP(JFLD+IOFF) = ISHFUV_G(ISTUV_G+JFLD-1)
        IOFF = IOFF+IF_UV_G
      ENDIF
      IF(LDIVGP) THEN
        IPTRGP(JFLD+IOFF) = IOFFD+ISHFUV_G(ISTUV_G+JFLD-1)
        IOFF = IOFF+IF_UV_G
      ENDIF
      IPTRGP(JFLD+IOFF) = IOFFU+ISHFUV_G(ISTUV_G+JFLD-1)
      IOFF = IOFF+IF_UV_G
      IPTRGP(JFLD+IOFF) = IOFFV+ISHFUV_G(ISTUV_G+JFLD-1)
      IOFF = IOFF+IF_UV_G+IF_SCALARS_G
      IF(LSCDERS) THEN
        IOFF  = IOFF+IF_SCALARS_G
      ENDIF
      IF(LUVDER) THEN
        IPTRGP(JFLD+IOFF) = IOFFUVD+ISHFUV_G(ISTUV_G+JFLD-1)
        IOFF = IOFF+IF_UV_G
        IPTRGP(JFLD+IOFF) = IOFFUVD+KF_UV_G+ISHFUV_G(ISTUV_G+JFLD-1)
      ENDIF
    ENDDO

    DO JFLD=1,IF_SCALARS_G
      IOFF = 2*IF_UV_G
      IF (LVORGP) IOFF = IOFF+IF_UV_G
      IF (LDIVGP) IOFF = IOFF+IF_UV_G
      IPTRGP(JFLD+IOFF) = IOFFSC+ISHFSC_G(ISTSC_G+JFLD-1)
      IOFF = IOFF+IF_SCALARS_G
      IF(LSCDERS) THEN
        IPTRGP(JFLD+IOFF) = IOFFSCNS+ISHFSC_G(ISTSC_G+JFLD-1)
        IOFF  = IOFF+IF_SCALARS_G
        IF(LUVDER) THEN
          IOFF  = IOFF+2*IF_UV_G
        ENDIF
        IPTRGP(JFLD+IOFF) = IOFFSCEW+ISHFSC_G(ISTSC_G+JFLD-1)
      ENDIF
    ENDDO
    DO JFLD=1,IF_UV
      IPTRSPUV(JFLD) = ISTUV+JFLD-1
    ENDDO
    DO JFLD=1,IF_SCALARS
      IPTRSPSC(JFLD) = ISTSC+JFLD-1
    ENDDO

    IF(IF_UV_G > 0 .AND. IF_SCALARS_G > 0) THEN
      CALL FTINV_CTLAD(PGP,IF_UV_G,IF_SCALARS_G,&
       & IF_UV,IF_SCALARS,IF_SCDERS,IF_GP,IF_FS,IF_OUT_LT,&
       & KVSETUV=IVSETUV(ISTUV_G:IENUV_G),&
       & KVSETSC=IVSETSC(ISTSC_G:IENSC_G),KPTRGP=IPTRGP)
    ELSEIF(IF_UV_G > 0) THEN
      CALL FTINV_CTLAD(PGP,IF_UV_G,IF_SCALARS_G,&
       & IF_UV,IF_SCALARS,IF_SCDERS,IF_GP,IF_FS,IF_OUT_LT,&
       & KVSETUV=IVSETUV(ISTUV_G:IENUV_G),&
       & KPTRGP=IPTRGP)
    ELSEIF(IF_SCALARS_G > 0) THEN
      CALL FTINV_CTLAD(PGP,IF_UV_G,IF_SCALARS_G,&
       & IF_UV,IF_SCALARS,IF_SCDERS,IF_GP,IF_FS,IF_OUT_LT,&
       & KVSETSC=IVSETSC(ISTSC_G:IENSC_G),KPTRGP=IPTRGP)
    ENDIF
    CALL LTINV_CTLAD(IF_OUT_LT,IF_UV,IF_SCALARS,IF_SCDERS, &
     & PSPVOR=PSPVOR,PSPDIV=PSPDIV,PSPSCALAR=PSPSCALAR,&
     & KFLDPTRUV=IPTRSPUV,KFLDPTRSC=IPTRSPSC)
    
  ENDDO

ELSE

  ! No splitting of fields, transform done in one go

  CALL FTINV_CTLAD(PGP,KF_UV_G,KF_SCALARS_G,&
   & KF_UV,KF_SCALARS,KF_SCDERS,KF_GP,KF_FS,KF_OUT_LT,&
   & KVSETUV,KVSETSC)

  CALL LTINV_CTLAD(KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS, &
   &PSPVOR=PSPVOR,PSPDIV=PSPDIV,PSPSCALAR=PSPSCALAR)

ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE INV_TRANS_CTLAD
END MODULE INV_TRANS_CTLAD_MOD
