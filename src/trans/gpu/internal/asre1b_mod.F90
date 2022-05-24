! (C) Copyright 2000- ECMWF.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE ASRE1B_MOD
CONTAINS
SUBROUTINE ASRE1B(KFIELD,PAOA,PSOA)

USE PARKIND_ECTRANS  ,ONLY : JPIM     ,JPRBT

USE TPM_DIM         ,ONLY : R, R_NDGNH, R_NDGL
USE TPM_TRANS       ,ONLY : FOUBUF_IN
USE TPM_GEOMETRY    ,ONLY : G, G_NDGLU
USE TPM_DISTR       ,ONLY : D,D_NUMP,D_MYMS,D_NSTAGT1B,D_NPROCL,D_NPNTGTB1
use tpm_gen, only: nout

!**** *ASRE1B* - Recombine antisymmetric and symmetric parts

!     Purpose.
!     --------
!        To recombine the antisymmetric and symmetric parts of the
!        Fourier arrays and update the correct parts of the state
!        variables.

!**   Interface.
!     ----------
!        *CALL* *ASRE1B(..)

!        Explicit arguments :
!        -------------------   KFIELD - number of fields (input-c)
!                              KM - zonal wavenumber(input-c)
!                              KMLOC - local version of KM (input-c)
!                              PAOA - antisymmetric part of Fourier
!                              fields for zonal wavenumber KM (input)
!                              PSOA - symmetric part of Fourier
!                              fields for zonal wavenumber KM (input)

!        Implicit arguments :  FOUBUF_IN - output buffer (output)
!        --------------------

!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From ASRE1B in IFS CY22R1

!     ------------------------------------------------------------------


IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KFIELD
INTEGER(KIND=JPIM) :: KM,KMLOC
REAL(KIND=JPRBT),   INTENT(IN)  :: PSOA(:,:,:)
REAL(KIND=JPRBT),   INTENT(IN)  :: PAOA(:,:,:)
!INTEGER(KIND=JPIM), INTENT(OUT) :: ISTAN(:,:)
!INTEGER(KIND=JPIM), INTENT(OUT) :: ISTAS(:,:)

!     LOCAL INTEGERS
INTEGER(KIND=JPIM) :: ISL, IGLS, JFLD, JGL ,IPROC, IPROCS, IDGNH, ISTAN, ISTAS

!     ------------------------------------------------------------------
 
!*       1.    RECOMBINATION  OF SYMMETRIC AND ANTSYMMETRIC PARTS.
!              ---------------------------------------------------

!$ACC DATA PRESENT(PAOA,PSOA,D_MYMS,D_NPROCL,D_NSTAGT1B,D_NPNTGTB1,G_NDGLU) PRESENT(FOUBUF_IN)
!$ACC PARALLEL LOOP DEFAULT(NONE) PRIVATE(KM,ISL,IPROC,ISTAN,IGLS,IPROCS,ISTAS)
DO KMLOC=1,D_NUMP
  KM = D_MYMS(KMLOC)
  ISL = MAX(R_NDGNH-G_NDGLU(KM)+1,1)
  DO JGL=ISL, R_NDGNH
!        if (JGL .ge. ISL) then
          IPROC = D_NPROCL(JGL)
          ISTAN = (D_NSTAGT1B(IPROC) + D_NPNTGTB1(KMLOC,JGL))*2*KFIELD
          IGLS = R_NDGL+1-JGL
          IPROCS = D_NPROCL(IGLS)
          ISTAS = (D_NSTAGT1B(IPROCS) + D_NPNTGTB1(KMLOC,IGLS))*2*KFIELD
          DO JFLD=1,2*KFIELD
            !write(iunit,*) 'xx ', KM, KFIELD, ISL , KMLOC, JFLD, ISTAN, ISTAS, IGLS, JGL, IPROC, PAOA(JFLD,JGL,KMLOC), PSOA(JFLD,JGL,KMLOC)
            !call flush(iunit)
            FOUBUF_IN(ISTAN+JFLD) = PAOA(JFLD,JGL,KMLOC)+PSOA(JFLD,JGL,KMLOC)
            FOUBUF_IN(ISTAS+JFLD) = PSOA(JFLD,JGL,KMLOC)-PAOA(JFLD,JGL,KMLOC)
          ENDDO
!        end if
   ENDDO
ENDDO
!$ACC END DATA

!     ------------------------------------------------------------------

END SUBROUTINE ASRE1B
END MODULE ASRE1B_MOD
