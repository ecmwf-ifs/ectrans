! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TPM_FFTW
!   Author.
!   -------
!     George Mozdzynski
!
!   Modifications.
!   -------------- 
!     Original      October 2014
!     R. El Khatib  01-Sep-2015 More subroutines for better modularity
!     R. El Khatib  08-Jun-2023 LALL_FFTW for better flexibility
!     W. Deconinck  17-Jun-2024 Replace legacy FFTW interface with the FFTW3 interface, add documentation, and improve clarity

USE, INTRINSIC :: ISO_C_BINDING

USE PARKIND1        ,ONLY : JPIM, JPRB, JPRD
USE MPL_MODULE      ,ONLY : MPL_MYRANK
USE YOMHOOK         ,ONLY : LHOOK, DR_HOOK, JPHOOK
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS

IMPLICIT NONE

SAVE

#ifdef __NEC__
! From NLC (NEC Numeric Library Collection)
#include "aslfftw3.f03"
#define FFTW_NO_SIMD 0
#else
#include "fftw3.f03"
#endif

! For now we use the LEGACY FFTW INTERFACE still, to be removed once we are sure no
! FFTW-pretending libraries are only implementing the legacy interface.
#define LEGACY_FFTW_INTERFACE 1

PRIVATE
PUBLIC INIT_PLANS_FFTW, DESTROY_PLANS_FFTW, FFTW_RESOL, TW, EXEC_FFTW, EXEC_EFFTW

!> @brief Cache state and execution settings for FFTW plans at one active resolution.
TYPE FFTW_TYPE
  !! Number of cached plans currently stored for each transform length `KN`.
  INTEGER(KIND=JPIM),ALLOCATABLE :: N_PLANS(:)

  !! Head nodes of the linked lists storing cached plans for each transform length.
  TYPE(FFTW_PLAN),POINTER :: FFTW_PLANS(:) => NULL()

  !! Largest transform length for which this cache has been initialised.
  INTEGER(KIND=JPIM) :: N_MAX=0

  !! N_MAX_PLANS: Maximum number of cached plans retained for any one transform length. The practical cache key is (KN, KTYPE, KLOT).
  !! EXEC_FFTW passes KN = KRLEN and KLOT = NBATCH, with NBATCH = 1 when LALL_FFTW=.FALSE. and NBATCH = KFIELDS when LALL_FFTW=.TRUE.
  !! Since KTYPE is only 1 or -1, eviction starts when one fixed KRLEN needs a fifth distinct (KTYPE, NBATCH) combination, when N_MAX_PLANS=4.
  INTEGER(KIND=JPIM) :: N_MAX_PLANS=4

  !! If true, execute all fields in one batched FFTW call; otherwise execute one field at a time.
  LOGICAL            :: LALL_FFTW=.FALSE.
END TYPE FFTW_TYPE

INTEGER(KIND=JPIM), PARAMETER :: NTYPE_C2R=1, NTYPE_R2C=-1
INTEGER(KIND=JPIM), PARAMETER :: TPM_FFTW_PLAN_FLAGS = FFTW_ESTIMATE + FFTW_NO_SIMD
INTEGER(KIND=JPIM), PARAMETER :: NPLAN_ID_UNINITIALISED = 123456
INTEGER(KIND=JPIM), PARAMETER :: NPLAN_ID_DESTROYED     = 999999

!> @brief Metadata for one cached FFTW plan together with its linked-list successor.
TYPE FFTW_PLAN
  !! Sentinel used to detect uninitialised or recycled plan records.
  INTEGER(KIND=JPIM) :: NPLAN_ID=NPLAN_ID_UNINITIALISED

  !! FFTW plan handle returned by the modern FFTW interface.
  TYPE(C_PTR)        :: NPLAN=C_NULL_PTR

  !! Batch size used when the plan was created.
  INTEGER(KIND=JPIM) :: NLOT

  !! Transform direction selector associated with the plan (`1` or `-1`).
  INTEGER(KIND=JPIM) :: NTYPE

  !! Pointer to the next cached plan for the same transform length.
  TYPE(FFTW_PLAN),POINTER :: NEXT_PLAN => NULL()
END TYPE FFTW_PLAN

!> @brief Array of per-resolution FFTW caches indexed by the active transform handle.
TYPE(FFTW_TYPE),ALLOCATABLE,TARGET :: FFTW_RESOL(:)
!> @brief Pointer to the FFTW cache associated with the currently active resolution.
TYPE(FFTW_TYPE),POINTER     :: TW

! Redefine JPCB as some FFTW implementations already define it in their fftw3.f03
#define JPCB TPM_FFTW_JPCB
INTEGER, PARAMETER :: JPCB = MERGE(C_DOUBLE_COMPLEX, C_FLOAT_COMPLEX, JPRB == JPRD)


INTERFACE TPM_FFTW_PLAN_MANY_DFT_C2R
#if LEGACY_FFTW_INTERFACE
  PROCEDURE TPM_FFTW_PLAN_MANY_DFT_C2R_LEGACY
  PROCEDURE TPM_FFTWF_PLAN_MANY_DFT_C2R_LEGACY
#else
  PROCEDURE FFTW_PLAN_MANY_DFT_C2R   ! defined in fftw3.f03
  PROCEDURE FFTWF_PLAN_MANY_DFT_C2R  ! defined in fftw3.f03
#endif
END INTERFACE TPM_FFTW_PLAN_MANY_DFT_C2R

INTERFACE TPM_FFTW_PLAN_MANY_DFT_R2C
#if LEGACY_FFTW_INTERFACE
  PROCEDURE TPM_FFTW_PLAN_MANY_DFT_R2C_LEGACY
  PROCEDURE TPM_FFTWF_PLAN_MANY_DFT_R2C_LEGACY
#else
  PROCEDURE FFTW_PLAN_MANY_DFT_R2C  ! defined in fftw3.f03
  PROCEDURE FFTWF_PLAN_MANY_DFT_R2C ! defined in fftw3.f03
#endif
END INTERFACE TPM_FFTW_PLAN_MANY_DFT_R2C

INTERFACE TPM_FFTW_EXECUTE_DFT_C2R
#if LEGACY_FFTW_INTERFACE
  PROCEDURE TPM_FFTW_EXECUTE_DFT_C2R_LEGACY
  PROCEDURE TPM_FFTWF_EXECUTE_DFT_C2R_LEGACY
#else
  PROCEDURE FFTW_EXECUTE_DFT_C2R    ! defined in fftw3.f03
  PROCEDURE FFTWF_EXECUTE_DFT_C2R   ! defined in fftw3.f03
#endif
  PROCEDURE TPM_FFTW_EXECUTE_DFT_C2R_RANK2
END INTERFACE TPM_FFTW_EXECUTE_DFT_C2R

INTERFACE TPM_FFTW_EXECUTE_DFT_R2C
#if LEGACY_FFTW_INTERFACE
  PROCEDURE TPM_FFTW_EXECUTE_DFT_R2C_LEGACY
  PROCEDURE TPM_FFTWF_EXECUTE_DFT_R2C_LEGACY
#else
  PROCEDURE FFTW_EXECUTE_DFT_R2C   ! defined in fftw3.f03
  PROCEDURE FFTWF_EXECUTE_DFT_R2C  ! defined in fftw3.f03
#endif
  PROCEDURE TPM_FFTW_EXECUTE_DFT_R2C_RANK2
END INTERFACE TPM_FFTW_EXECUTE_DFT_R2C

! ------------------------------------------------------------------
CONTAINS
! ------------------------------------------------------------------

FUNCTION TPM_FFTW_ALLOC_COMPLEX(N) RESULT(RES)
!> @brief Allocate a complex work buffer through the precision-specific FFTW backend.
!!
!! @param[in] N Number of complex elements to allocate.
IMPLICIT NONE
INTEGER(KIND=C_SIZE_T),INTENT(IN) :: N
TYPE(C_PTR) :: RES
IF (JPRB == JPRD) THEN
  RES=FFTW_ALLOC_COMPLEX(N)
ELSE
  RES=FFTWF_ALLOC_COMPLEX(N)
END IF
END FUNCTION TPM_FFTW_ALLOC_COMPLEX


SUBROUTINE TPM_FFTW_FREE(PTR)
!> @brief Release a complex work buffer allocated through the FFTW precision wrapper.
!!
!! @param[in] PTR C pointer to a work buffer allocated by `TPM_FFTW_ALLOC_COMPLEX`.
IMPLICIT NONE
TYPE(C_PTR),INTENT(IN) :: PTR
IF (JPRB == JPRD) THEN
  CALL FFTW_FREE(PTR)
ELSE
  CALL FFTWF_FREE(PTR)
END IF
END SUBROUTINE TPM_FFTW_FREE


SUBROUTINE TPM_FFTW_DESTROY_PLAN(PLAN)
!> @brief Destroy an FFTW plan using the precision-specific backend.
!!
!! @param[in] PLAN FFTW plan handle to destroy.
IMPLICIT NONE
TYPE(C_PTR),INTENT(IN) :: PLAN
#if LEGACY_FFTW_INTERFACE == 1
INTEGER(KIND=C_INTPTR_T) :: PLAN_LEGACY
PLAN_LEGACY=TPM_FFTW_LEGACY_PLAN_FROM_CPTR(PLAN)
IF (JPRB == JPRD) THEN
  CALL DFFTW_DESTROY_PLAN(PLAN_LEGACY)
ELSE
  CALL SFFTW_DESTROY_PLAN(PLAN_LEGACY)
END IF
#else
IF (JPRB == JPRD) THEN
  CALL FFTW_DESTROY_PLAN(PLAN)
ELSE
  CALL FFTWF_DESTROY_PLAN(PLAN)
END IF
#endif
END SUBROUTINE TPM_FFTW_DESTROY_PLAN


SUBROUTINE TPM_FFTW_EXECUTE_DFT_C2R_RANK2(PLAN,IN,OUT)
!> @brief Execute a complex-to-real FFTW plan on rank-2 arrays via flattened views.
!!
!! @param[in] PLAN FFTW plan handle to execute.
!! @param[inout] IN Rank-2 complex work array passed to FFTW as a contiguous rank-1 view.
!! @param[inout] OUT Rank-2 real work array passed to FFTW as a contiguous rank-1 view.
IMPLICIT NONE
TYPE(C_PTR),INTENT(IN) :: PLAN
COMPLEX(KIND=JPCB),INTENT(INOUT),CONTIGUOUS,TARGET :: IN(:,:)
REAL(KIND=JPRB),INTENT(INOUT),CONTIGUOUS,TARGET :: OUT(:,:)
COMPLEX(KIND=JPCB),POINTER :: IN_FLAT(:)
REAL(KIND=JPRB),POINTER :: OUT_FLAT(:)
IN_FLAT(1:SIZE(IN)) => IN
OUT_FLAT(1:SIZE(OUT)) => OUT
CALL TPM_FFTW_EXECUTE_DFT_C2R(PLAN,IN_FLAT,OUT_FLAT)
END SUBROUTINE TPM_FFTW_EXECUTE_DFT_C2R_RANK2


SUBROUTINE TPM_FFTW_EXECUTE_DFT_R2C_RANK2(PLAN,IN,OUT)
!> @brief Execute a real-to-complex FFTW plan on rank-2 arrays via flattened views.
!!
!! @param[in] PLAN FFTW plan handle to execute.
!! @param[inout] IN Rank-2 real work array passed to FFTW as a contiguous rank-1 view.
!! @param[inout] OUT Rank-2 complex work array passed to FFTW as a contiguous rank-1 view.
IMPLICIT NONE
TYPE(C_PTR),INTENT(IN) :: PLAN
REAL(KIND=JPRB),INTENT(INOUT),CONTIGUOUS,TARGET :: IN(:,:)
COMPLEX(KIND=JPCB),INTENT(INOUT),CONTIGUOUS,TARGET :: OUT(:,:)
REAL(KIND=JPRB),POINTER :: IN_FLAT(:)
COMPLEX(KIND=JPCB),POINTER :: OUT_FLAT(:)
IN_FLAT(1:SIZE(IN)) => IN
OUT_FLAT(1:SIZE(OUT)) => OUT
CALL TPM_FFTW_EXECUTE_DFT_R2C(PLAN,IN_FLAT,OUT_FLAT)
END SUBROUTINE TPM_FFTW_EXECUTE_DFT_R2C_RANK2


SUBROUTINE INIT_PLANS_FFTW(KDLON)
!> @brief Allocate plan bookkeeping for the active resolution.
!!
!! @param[in] KDLON Upper bound on the transform lengths indexed in the plan cache.
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN) :: KDLON

TW%N_MAX=KDLON
ALLOCATE(TW%FFTW_PLANS(TW%N_MAX))
ALLOCATE(TW%N_PLANS(TW%N_MAX))
TW%N_PLANS(:)=0
RETURN  
END SUBROUTINE INIT_PLANS_FFTW


SUBROUTINE CREATE_PLAN_FFTW(KPLAN,KTYPE,KN,KLOT)
!> @brief Reuse or create a cached FFTW plan for a given transform shape.
!!
!! @param[out] KPLAN FFTW plan handle matching the requested transform layout.
!! @param[in] KTYPE Transform direction selector: `1` for complex-to-real and `-1` for real-to-complex.
!! @param[in] KN Real transform length used as the one-dimensional FFT extent.
!! @param[in] KLOT Number of transforms to execute in the batched FFTW plan.
!! @note Access to the shared plan cache is serialized by the named OpenMP critical region `FFTW_CREATE`.
!!   This includes plan lookup, eviction, and creation, so heavy multi-threaded contention here can limit
!!   scalability unless the needed plans are created early and then reused.
IMPLICIT NONE
TYPE(C_PTR),INTENT(OUT) :: KPLAN
INTEGER(KIND=JPIM),INTENT(IN) :: KTYPE,KN,KLOT

TYPE(C_PTR) :: IPLAN
INTEGER(KIND=JPIM) :: IRANK, ISTRIDE
INTEGER(KIND=JPIM) :: JL
INTEGER(KIND=JPIM) :: IRDIST,ICDIST,IN(1),IEMBED(1)
INTEGER(KIND=JPIM) :: CEMBED(1)
REAL(KIND=JPRB), POINTER :: ZDUM(:)
COMPLEX(KIND=JPCB), POINTER :: CDUM(:)
TYPE(C_PTR) :: ZDUMP
LOGICAL :: LLFOUND
LOGICAL, PARAMETER :: LLRESTRICT_PLANS=.TRUE.
TYPE(FFTW_PLAN),POINTER :: CURR_FFTW_PLAN,START_FFTW_PLAN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE, ZHOOK_HANDLE2
IF (LHOOK) CALL DR_HOOK('CREATE_PLAN_FFTW',0,ZHOOK_HANDLE)

IF( KN > TW%N_MAX )THEN
  CALL ABORT_TRANS('CREATE_PLAN_FFTW: KN > N_MAX THAT WAS INITIALISED IN INIT_PLANS_FFTW')
ENDIF

IRANK=1
ISTRIDE=1
IN(1)=KN
IEMBED(1)=IN(1)
ICDIST=KN/2+1
CEMBED(1)=ICDIST
IRDIST=ICDIST*2

!$OMP CRITICAL (FFTW_CREATE)
LLFOUND=.FALSE.
IF( TW%FFTW_PLANS(KN)%NPLAN_ID /= NPLAN_ID_UNINITIALISED )THEN
  WRITE(*,'("CREATE_PLAN_FFTW.1: PLAN_ID=",I10)')TW%FFTW_PLANS(KN)%NPLAN_ID
  CALL ABORT_TRANS('CREATE_PLAN_FFTW.1: NPLAN_ID /= NPLAN_ID_UNINITIALISED')
ENDIF
CURR_FFTW_PLAN=>TW%FFTW_PLANS(KN)
IF( CURR_FFTW_PLAN%NPLAN_ID /= NPLAN_ID_UNINITIALISED )THEN
  WRITE(*,'("CREATE_PLAN_FFTW.2: PLAN_ID=",I10)')CURR_FFTW_PLAN%NPLAN_ID
  CALL ABORT_TRANS('CREATE_PLAN_FFTW.2: NPLAN_ID /= NPLAN_ID_UNINITIALISED')
ENDIF
! search for plan in existing plans
DO JL=1,TW%N_PLANS(KN)
  IF( KLOT == CURR_FFTW_PLAN%NLOT .AND. KTYPE == CURR_FFTW_PLAN%NTYPE )THEN
    LLFOUND=.TRUE.
    IPLAN=CURR_FFTW_PLAN%NPLAN
    EXIT
  ELSEIF( JL /= TW%N_PLANS(KN) )THEN
    CURR_FFTW_PLAN=>CURR_FFTW_PLAN%NEXT_PLAN
    IF( CURR_FFTW_PLAN%NPLAN_ID /= NPLAN_ID_UNINITIALISED )THEN
      WRITE(*,'("CREATE_PLAN_FFTW.3: PLAN_ID=",I10)')CURR_FFTW_PLAN%NPLAN_ID
      CALL ABORT_TRANS('CREATE_PLAN_FFTW.3: NPLAN_ID /= NPLAN_ID_UNINITIALISED')
    ENDIF
  ENDIF
ENDDO
IF( .NOT.LLFOUND )THEN
  IF( LLRESTRICT_PLANS )THEN
    IF( TW%N_PLANS(KN) == TW%N_MAX_PLANS )THEN
      ! destroy the plan at the start of the list
      CALL TPM_FFTW_DESTROY_PLAN(TW%FFTW_PLANS(KN)%NPLAN)
      TW%FFTW_PLANS(KN)%NPLAN_ID=NPLAN_ID_DESTROYED
       ! mark the plan record as destroyed but keep it in the list to be recycled by the next plan creation with the same KN,
       ! which avoids costly deallocation and reallocation of the FFTW_PLAN record and preserves the linked list structure
       ! without needing to update any pointers
      START_FFTW_PLAN=>TW%FFTW_PLANS(KN)
      TW%FFTW_PLANS(KN)=TW%FFTW_PLANS(KN)%NEXT_PLAN
      ! DEALLOCATE(START_FFTW_PLAN)
      TW%N_PLANS(KN)=TW%N_PLANS(KN)-1
    ENDIF
  ENDIF
  ZDUMP=TPM_FFTW_ALLOC_COMPLEX(INT(1,C_SIZE_T))
  CALL C_F_POINTER(ZDUMP,ZDUM,[2])
  CALL C_F_POINTER(ZDUMP,CDUM,[1])
  IF( KTYPE==NTYPE_C2R )THEN
    IF (LHOOK) CALL DR_HOOK('FFTW_PLAN_MANY_DFT_C2R',0,ZHOOK_HANDLE2)
    IPLAN=TPM_FFTW_PLAN_MANY_DFT_C2R(IRANK,IN,KLOT,CDUM,CEMBED,ISTRIDE,ICDIST,&
         & ZDUM,IEMBED,ISTRIDE,IRDIST,TPM_FFTW_PLAN_FLAGS)
    IF (LHOOK) CALL DR_HOOK('FFTW_PLAN_MANY_DFT_C2R',1,ZHOOK_HANDLE2)
  ELSEIF( KTYPE==NTYPE_R2C )THEN
    IF (LHOOK) CALL DR_HOOK('FFTW_PLAN_MANY_DFT_R2C',0,ZHOOK_HANDLE2)
    IPLAN=TPM_FFTW_PLAN_MANY_DFT_R2C(IRANK,IN,KLOT,ZDUM,IEMBED,ISTRIDE,IRDIST,&
         & CDUM,CEMBED,ISTRIDE,ICDIST,TPM_FFTW_PLAN_FLAGS)
    IF (LHOOK) CALL DR_HOOK('FFTW_PLAN_MANY_DFT_R2C',1,ZHOOK_HANDLE2)
  ELSE
    CALL ABORT_TRANS('FFTW_PLAN: INVALID KTYPE')
  ENDIF
  CALL TPM_FFTW_FREE(ZDUMP)
  KPLAN=IPLAN
  TW%N_PLANS(KN)=TW%N_PLANS(KN)+1
  IF( TW%N_PLANS(KN) /= 1 )THEN
    ALLOCATE(CURR_FFTW_PLAN%NEXT_PLAN)
    CURR_FFTW_PLAN=>CURR_FFTW_PLAN%NEXT_PLAN
  ENDIF
  IF( CURR_FFTW_PLAN%NPLAN_ID /= NPLAN_ID_UNINITIALISED )THEN
    WRITE(*,'("CREATE_PLAN_FFTW.4: PLAN_ID=",I10)')CURR_FFTW_PLAN%NPLAN_ID
    CALL ABORT_TRANS('CREATE_PLAN_FFTW.4: NPLAN_ID /= NPLAN_ID_UNINITIALISED')
  ENDIF
  CURR_FFTW_PLAN%NPLAN=IPLAN
  CURR_FFTW_PLAN%NLOT=KLOT
  CURR_FFTW_PLAN%NTYPE=KTYPE
  CURR_FFTW_PLAN%NEXT_PLAN=>NULL()
ELSE
  KPLAN=IPLAN
ENDIF
!$OMP END CRITICAL (FFTW_CREATE)

IF (LHOOK) CALL DR_HOOK('CREATE_PLAN_FFTW',1,ZHOOK_HANDLE)
RETURN
END SUBROUTINE CREATE_PLAN_FFTW


SUBROUTINE DESTROY_PLAN_FFTW(KPLAN)
!> @brief Destroy one cached FFTW plan inside the module-level OpenMP critical region.
!!
!! @param[in] KPLAN FFTW plan handle to destroy.
!! @note Plan destruction is serialized by the named OpenMP critical region `FFTW_DESTROY`.
!!   This is usually a teardown path, but repeated destruction from many threads will not run concurrently.
IMPLICIT NONE
TYPE(C_PTR),INTENT(IN) :: KPLAN
!$OMP CRITICAL (FFTW_DESTROY)
CALL TPM_FFTW_DESTROY_PLAN(KPLAN)
!$OMP END CRITICAL (FFTW_DESTROY)
RETURN
END SUBROUTINE DESTROY_PLAN_FFTW


SUBROUTINE DESTROY_PLANS_FFTW
!> @brief Destroy all cached FFTW plans and release the resolution-local plan tables.
!!
!! @note Plans are destroyed one by one through `DESTROY_PLAN_FFTW`, so teardown repeatedly enters the
!!   `FFTW_DESTROY` critical region and is serialized across threads.
IMPLICIT NONE
INTEGER(KIND=JPIM) :: JL, JN
TYPE(FFTW_PLAN),POINTER :: CURR_FFTW_PLAN, NEXT_FFTW_PLAN
DO JN=1,TW%N_MAX
  CURR_FFTW_PLAN=>TW%FFTW_PLANS(JN)
  DO JL=1,TW%N_PLANS(JN)
    CALL DESTROY_PLAN_FFTW(CURR_FFTW_PLAN%NPLAN)
    NEXT_FFTW_PLAN=>CURR_FFTW_PLAN%NEXT_PLAN
    IF( JL /= 1 ) THEN
      DEALLOCATE( CURR_FFTW_PLAN )
    ENDIF
    CURR_FFTW_PLAN => NEXT_FFTW_PLAN
  ENDDO
ENDDO
IF( ASSOCIATED(TW) ) THEN
  IF( ASSOCIATED(TW%FFTW_PLANS) ) DEALLOCATE(TW%FFTW_PLANS)
  IF( ALLOCATED(TW%N_PLANS) )     DEALLOCATE(TW%N_PLANS)
  TW%N_MAX=0
ENDIF
RETURN
END SUBROUTINE DESTROY_PLANS_FFTW


SUBROUTINE COPY_PREEL_TO_ZFFT(KLEN,KOFF,LD_TRANSPOSED,PREEL,ZFFT)
!> @brief Copy a batched `PREEL` slab into the rank-2 FFT work array.
!!
!! @param[in] KLEN Number of values copied per field.
!! @param[in] KOFF One-based offset of the transform segment inside `PREEL`.
!! @param[in] LD_TRANSPOSED Selects `PREEL(point,field)` when true and `PREEL(field,point)` when false.
!! @param[in] PREEL Source work array.
!! @param[inout] ZFFT Destination rank-2 FFT work buffer.
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN) :: KLEN
INTEGER(KIND=JPIM),INTENT(IN) :: KOFF
LOGICAL,INTENT(IN) :: LD_TRANSPOSED
REAL(KIND=JPRB),INTENT(IN) :: PREEL(:,:)
REAL(KIND=JPRB),INTENT(INOUT) :: ZFFT(:,:)

INTEGER(KIND=JPIM) :: JJ,JF

IF( LD_TRANSPOSED )THEN
  DO JF=1,SIZE(ZFFT,2)
    DO JJ=1,KLEN
      ZFFT(JJ,JF)=PREEL(KOFF+JJ-1,JF)
    ENDDO
  ENDDO
ELSE
  DO JF=1,SIZE(ZFFT,2)
    DO JJ=1,KLEN
      ZFFT(JJ,JF)=PREEL(JF,KOFF+JJ-1)
    ENDDO
  ENDDO
ENDIF
END SUBROUTINE COPY_PREEL_TO_ZFFT

SUBROUTINE COPY_ZFFT_TO_PREEL(KLEN,KOFF,LD_TRANSPOSED,ZFFT,PREEL,PSCALE)
!> @brief Copy a rank-2 FFT work array back into `PREEL`, optionally applying a scale factor.
!!
!! @param[in] KLEN Number of values copied per field.
!! @param[in] KOFF One-based offset of the transform segment inside `PREEL`.
!! @param[in] LD_TRANSPOSED Selects `PREEL(point,field)` when true and `PREEL(field,point)` when false.
!! @param[in] ZFFT Source rank-2 FFT work buffer.
!! @param[inout] PREEL Destination work array.
!! @param[in] PSCALE Optional scale factor applied to each copied value when present.
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN) :: KLEN
INTEGER(KIND=JPIM),INTENT(IN) :: KOFF
LOGICAL,INTENT(IN) :: LD_TRANSPOSED
REAL(KIND=JPRB),INTENT(IN) :: ZFFT(:,:)
REAL(KIND=JPRB),INTENT(INOUT) :: PREEL(:,:)
REAL(KIND=JPRB), OPTIONAL,INTENT(IN) :: PSCALE

INTEGER(KIND=JPIM) :: JJ,JF

IF( LD_TRANSPOSED )THEN
  IF (PRESENT(PSCALE)) THEN
     DO JF=1,SIZE(ZFFT,2)
       DO JJ=1,KLEN
         PREEL(KOFF+JJ-1,JF)=ZFFT(JJ,JF)*PSCALE
       ENDDO
     ENDDO
  ELSE
    DO JF=1,SIZE(ZFFT,2)
      DO JJ=1,KLEN
        PREEL(KOFF+JJ-1,JF)=ZFFT(JJ,JF)
      ENDDO
    ENDDO
  ENDIF
ELSE
  IF (PRESENT(PSCALE)) THEN
     DO JF=1,SIZE(ZFFT,2)
       DO JJ=1,KLEN
         PREEL(JF,KOFF+JJ-1)=ZFFT(JJ,JF)*PSCALE
       ENDDO
     ENDDO
  ELSE
    DO JF=1,SIZE(ZFFT,2)
      DO JJ=1,KLEN
        PREEL(JF,KOFF+JJ-1)=ZFFT(JJ,JF)
      ENDDO
    ENDDO
  ENDIF
ENDIF
END SUBROUTINE COPY_ZFFT_TO_PREEL


SUBROUTINE COPY_PREEL_JF_TO_ZFFT_1(KLEN,KOFF,LD_TRANSPOSED,PREEL,KFIELD,ZFFT1)
!> @brief Copy one field from `PREEL` into the rank-1 FFT work array.
!!
!! @param[in] KLEN Number of values copied for the selected field.
!! @param[in] KOFF One-based offset of the transform segment inside `PREEL`.
!! @param[in] LD_TRANSPOSED Selects `PREEL(point,field)` when true and `PREEL(field,point)` when false.
!! @param[in] PREEL Source work array.
!! @param[in] KFIELD Field index copied into the work array.
!! @param[inout] ZFFT1 Destination rank-1 FFT work buffer.
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN) :: KLEN
INTEGER(KIND=JPIM),INTENT(IN) :: KOFF
LOGICAL,INTENT(IN) :: LD_TRANSPOSED
REAL(KIND=JPRB),INTENT(IN) :: PREEL(:,:)
INTEGER(KIND=JPIM),INTENT(IN) :: KFIELD
REAL(KIND=JPRB),INTENT(INOUT) :: ZFFT1(:)

INTEGER(KIND=JPIM) :: JJ

IF( LD_TRANSPOSED )THEN
  DO JJ=1,KLEN
    ZFFT1(JJ)=PREEL(KOFF+JJ-1,KFIELD)
  ENDDO
ELSE
  DO JJ=1,KLEN
    ZFFT1(JJ)=PREEL(KFIELD,KOFF+JJ-1)
  ENDDO
ENDIF
END SUBROUTINE COPY_PREEL_JF_TO_ZFFT_1


SUBROUTINE COPY_ZFFT_1_TO_PREEL_JF(KLEN,KOFF,LD_TRANSPOSED,ZFFT1,PREEL,KFIELD,PSCALE)
!> @brief Copy one FFT work vector back into `PREEL`, optionally applying a scale factor.
!!
!! @param[in] KLEN Number of values copied for the selected field.
!! @param[in] KOFF One-based offset of the transform segment inside `PREEL`.
!! @param[in] LD_TRANSPOSED Selects `PREEL(point,field)` when true and `PREEL(field,point)` when false.
!! @param[in] ZFFT1 Source rank-1 FFT work buffer.
!! @param[inout] PREEL Destination work array.
!! @param[in] KFIELD Field index updated in `PREEL`.
!! @param[in] PSCALE Optional scale factor applied to each copied value when present.
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN) :: KLEN
INTEGER(KIND=JPIM),INTENT(IN) :: KOFF
LOGICAL,INTENT(IN) :: LD_TRANSPOSED
REAL(KIND=JPRB),INTENT(IN) :: ZFFT1(:)
REAL(KIND=JPRB),INTENT(INOUT) :: PREEL(:,:)
INTEGER(KIND=JPIM),INTENT(IN) :: KFIELD
REAL(KIND=JPRB), OPTIONAL,INTENT(IN) :: PSCALE

INTEGER(KIND=JPIM) :: JJ

IF( LD_TRANSPOSED )THEN
  IF (PRESENT(PSCALE)) THEN
    DO JJ=1,KLEN
      PREEL(KOFF+JJ-1,KFIELD)=ZFFT1(JJ)*PSCALE
    ENDDO
  ELSE
    DO JJ=1,KLEN
      PREEL(KOFF+JJ-1,KFIELD)=ZFFT1(JJ)
    ENDDO
  ENDIF
ELSE
  IF (PRESENT(PSCALE)) THEN
    DO JJ=1,KLEN
      PREEL(KFIELD,KOFF+JJ-1)=ZFFT1(JJ)*PSCALE
    ENDDO
  ELSE
    DO JJ=1,KLEN
      PREEL(KFIELD,KOFF+JJ-1)=ZFFT1(JJ)
    ENDDO
  ENDIF
ENDIF
END SUBROUTINE COPY_ZFFT_1_TO_PREEL_JF


SUBROUTINE EXEC_FFTW_IMPL(CDNAME,KTYPE,KRLEN,KCLEN,KOFF,KFIELDS,LD_ALL,PREEL,LD_TRANSPOSED)
!> @brief Execute the shared FFTW path for both PREEL memory layouts and batching modes.
!!
!! @param[in] CDNAME Name passed to `DR_HOOK` and used in runtime error messages.
!! @param[in] KTYPE Transform direction selector: `1` for complex-to-real and `-1` for real-to-complex.
!! @param[in] KRLEN Number of real values in each transform.
!! @param[in] KCLEN Storage length of the packed FFT segment in real-valued form.
!! @param[in] KOFF One-based offset of the transform segment inside `PREEL`.
!! @param[in] KFIELDS Number of fields to transform.
!! @param[in] LD_ALL If true, execute all fields in one batched FFTW call; otherwise transform one field at a time.
!! @param[inout] PREEL Work array holding the transform input on entry and the transformed output on return.
!! @param[in] LD_TRANSPOSED Selects the `PREEL(point,field)` layout when true and `PREEL(field,point)` when false.
!! @note The internal call to `CREATE_PLAN_FFTW` accesses the shared FFTW plan cache inside the `FFTW_CREATE`
!!   critical region. Reusing warmed-up plans avoids repeated planning work, but cache access itself remains
!!   serialized across threads.
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN) :: CDNAME
INTEGER(KIND=JPIM),INTENT(IN) :: KTYPE
INTEGER(KIND=JPIM),INTENT(IN) :: KRLEN
INTEGER(KIND=JPIM),INTENT(IN) :: KCLEN
INTEGER(KIND=JPIM),INTENT(IN) :: KOFF
INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS
LOGICAL,INTENT(IN) :: LD_ALL
REAL(KIND=JPRB),INTENT(INOUT) :: PREEL(:,:)
LOGICAL,INTENT(IN) :: LD_TRANSPOSED

REAL(KIND=JPRB), POINTER :: ZFFT(:,:)
COMPLEX(KIND=JPCB), POINTER :: CFFT(:,:)
TYPE(C_PTR) :: ZFFTP
TYPE(C_PTR) :: IPLAN
INTEGER(KIND=JPIM) :: JF
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE, ZHOOK_HANDLE2
INTEGER(KIND=JPIM) :: NBATCH
REAL(KIND=JPRB) :: ZSCALE

IF (LHOOK) CALL DR_HOOK(CDNAME,0,ZHOOK_HANDLE)

IF ( (KTYPE /= NTYPE_R2C) .AND. (KTYPE /= NTYPE_C2R) ) THEN
  CALL ABORT_TRANS('TPM_FFTW:'//TRIM(CDNAME)//' : WRONG VALUE KTYPE')
ENDIF

NBATCH = MERGE(KFIELDS,1,LD_ALL)
CALL CREATE_PLAN_FFTW(IPLAN, KTYPE, KRLEN, NBATCH)

! Allocate a work array large enough to hold the complex FFT data for the entire batch of fields.
! The size is KCLEN/2 complex numbers per field, and there are NBATCH fields in the batch.
! This size is also sufficient for the real data
ZFFTP=TPM_FFTW_ALLOC_COMPLEX(INT(KCLEN/2*NBATCH,C_SIZE_T))

! It is chosen to alias the complex and real views of the work array through C_F_POINTER,
! so that FFTW performs in-place transforms.
! It should be investigated if out-of-place transforms with separate work arrays for real and complex data can improve performance
! by allowing more flexible memory access patterns in FFTW.
CALL C_F_POINTER(ZFFTP,ZFFT,[KCLEN,  NBATCH])
CALL C_F_POINTER(ZFFTP,CFFT,[KCLEN/2,NBATCH])

IF( LD_ALL ) THEN
  ! All fields are transformed together in a single FFTW call, so the work array is laid out as ZFFT(transform_point,field).
  IF (KTYPE==NTYPE_C2R) THEN
    CALL COPY_PREEL_TO_ZFFT(KCLEN,KOFF,LD_TRANSPOSED,PREEL,ZFFT)
    IF (LHOOK) CALL DR_HOOK('FFTW_EXECUTE_DFT_C2R',0,ZHOOK_HANDLE2)
    CALL TPM_FFTW_EXECUTE_DFT_C2R(IPLAN,CFFT,ZFFT)
    IF (LHOOK) CALL DR_HOOK('FFTW_EXECUTE_DFT_C2R',1,ZHOOK_HANDLE2)
    CALL COPY_ZFFT_TO_PREEL(KRLEN,KOFF,LD_TRANSPOSED,ZFFT,PREEL)
  ELSE
    CALL COPY_PREEL_TO_ZFFT(KRLEN,KOFF,LD_TRANSPOSED,PREEL,ZFFT)
    IF (LHOOK) CALL DR_HOOK('FFTW_EXECUTE_DFT_R2C',0,ZHOOK_HANDLE2)
    CALL TPM_FFTW_EXECUTE_DFT_R2C(IPLAN,ZFFT,CFFT)
    IF (LHOOK) CALL DR_HOOK('FFTW_EXECUTE_DFT_R2C',1,ZHOOK_HANDLE2)
    ! Real-to-complex transforms require scaling by 1/KRLEN, which can be applied in the copy back to PREEL to save cycles in the FFT work array.
    ZSCALE = 1.0_JPRB/REAL(KRLEN,JPRB)
    CALL COPY_ZFFT_TO_PREEL(KCLEN,KOFF,LD_TRANSPOSED,ZFFT,PREEL,PSCALE=ZSCALE)
  ENDIF
ELSE
  ! All fields are transformed separately in a loop over `JF`, so the work array is laid out as ZFFT(transform_point,1) and only one field is copied in and out at a time.
  IF (KTYPE==NTYPE_C2R) THEN
    DO JF=1,KFIELDS
      CALL COPY_PREEL_JF_TO_ZFFT_1(KCLEN,KOFF,LD_TRANSPOSED,PREEL,JF,ZFFT(:,1))
      IF (LHOOK) CALL DR_HOOK('FFTW_EXECUTE_DFT_C2R',0,ZHOOK_HANDLE2)
      CALL TPM_FFTW_EXECUTE_DFT_C2R(IPLAN,CFFT(:,1),ZFFT(:,1))
      IF (LHOOK) CALL DR_HOOK('FFTW_EXECUTE_DFT_C2R',1,ZHOOK_HANDLE2)
      CALL COPY_ZFFT_1_TO_PREEL_JF(KRLEN,KOFF,LD_TRANSPOSED,ZFFT(:,1),PREEL,JF)
    ENDDO
  ELSE
    ! Real-to-complex transforms require scaling by 1/KRLEN, which can be applied in the copy back to PREEL to save cycles in the FFT work array.
    ZSCALE = 1.0_JPRB/REAL(KRLEN,JPRB)
    DO JF=1,KFIELDS
      CALL COPY_PREEL_JF_TO_ZFFT_1(KRLEN,KOFF,LD_TRANSPOSED,PREEL,JF,ZFFT(:,1))
      IF (LHOOK) CALL DR_HOOK('FFTW_EXECUTE_DFT_R2C',0,ZHOOK_HANDLE2)
      CALL TPM_FFTW_EXECUTE_DFT_R2C(IPLAN,ZFFT(:,1),CFFT(:,1))
      IF (LHOOK) CALL DR_HOOK('FFTW_EXECUTE_DFT_R2C',1,ZHOOK_HANDLE2)
      CALL COPY_ZFFT_1_TO_PREEL_JF(KCLEN,KOFF,LD_TRANSPOSED,ZFFT(:,1),PREEL,JF,PSCALE=ZSCALE)
    ENDDO
  ENDIF
ENDIF
CALL TPM_FFTW_FREE(ZFFTP)

IF (LHOOK) CALL DR_HOOK(CDNAME,1,ZHOOK_HANDLE)
END SUBROUTINE EXEC_FFTW_IMPL


SUBROUTINE EXEC_FFTW(KTYPE,KRLEN,KCLEN,KOFF,KFIELDS,LD_ALL,PREEL)
!> @brief Execute FFTW transforms for arrays stored as `PREEL(field,point)`.
!!
!! @param[in] KTYPE Transform direction selector: `1` for complex-to-real and `-1` for real-to-complex.
!! @param[in] KRLEN Number of real values in each transform.
!! @param[in] KCLEN Storage length of the packed FFT segment in real-valued form.
!! @param[in] KOFF One-based offset of the transform segment inside `PREEL`.
!! @param[in] KFIELDS Number of fields to transform.
!! @param[in] LD_ALL If true, execute all fields in one batched FFTW call; otherwise transform one field at a time.
!! @param[inout] PREEL Work array laid out as `PREEL(field,point)` containing the input on entry and the output on return.
!! @note This routine shares the module-wide FFTW plan cache, so calls may serialize briefly in `CREATE_PLAN_FFTW`
!!   when accessing the `FFTW_CREATE` critical region.
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)   :: KTYPE
INTEGER(KIND=JPIM),INTENT(IN)   :: KRLEN
INTEGER(KIND=JPIM),INTENT(IN)   :: KCLEN
INTEGER(KIND=JPIM),INTENT(IN)   :: KOFF
INTEGER(KIND=JPIM),INTENT(IN)   :: KFIELDS
LOGICAL           ,INTENT(IN)   :: LD_ALL
REAL(KIND=JPRB), INTENT(INOUT)  :: PREEL(:,:)

CALL EXEC_FFTW_IMPL('EXEC_FFTW',KTYPE,KRLEN,KCLEN,KOFF,KFIELDS,LD_ALL,PREEL,LD_TRANSPOSED=.FALSE.)
END SUBROUTINE EXEC_FFTW

SUBROUTINE EXEC_EFFTW(KTYPE,KRLEN,KCLEN,KOFF,KFIELDS,LD_ALL,PREEL)
!> @brief Execute FFTW transforms for arrays stored as `PREEL(point,field)`.
!!
!! @param[in] KTYPE Transform direction selector: `1` for complex-to-real and `-1` for real-to-complex.
!! @param[in] KRLEN Number of real values in each transform.
!! @param[in] KCLEN Storage length of the packed FFT segment in real-valued form.
!! @param[in] KOFF One-based offset of the transform segment inside `PREEL`.
!! @param[in] KFIELDS Number of fields to transform.
!! @param[in] LD_ALL If true, execute all fields in one batched FFTW call; otherwise transform one field at a time.
!! @param[inout] PREEL Work array laid out as `PREEL(point,field)` containing the input on entry and the output on return.
!! @note This routine shares the module-wide FFTW plan cache, so calls may serialize briefly in `CREATE_PLAN_FFTW`
!!   when accessing the `FFTW_CREATE` critical region. Reusing warmed-up plans avoids repeated planning work,
!!   but cache access itself remains serialized across threads.
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)   :: KTYPE
INTEGER(KIND=JPIM),INTENT(IN)   :: KRLEN
INTEGER(KIND=JPIM),INTENT(IN)   :: KCLEN
INTEGER(KIND=JPIM),INTENT(IN)   :: KOFF
INTEGER(KIND=JPIM),INTENT(IN)   :: KFIELDS
LOGICAL           ,INTENT(IN)   :: LD_ALL
REAL(KIND=JPRB), INTENT(INOUT)  :: PREEL(:,:)

CALL EXEC_FFTW_IMPL('EXEC_EFFTW',KTYPE,KRLEN,KCLEN,KOFF,KFIELDS,LD_ALL,PREEL,LD_TRANSPOSED=.TRUE.)
END SUBROUTINE EXEC_EFFTW

! -----------------------------------------------------------------------------
! Following routines are for legacy FFTW interface support and can be removed
! once the legacy interface is no longer needed.
! -----------------------------------------------------------------------------

#if LEGACY_FFTW_INTERFACE
FUNCTION TPM_FFTW_CPTR_FROM_LEGACY_PLAN(PLAN_LEGACY) RESULT(PLAN)
IMPLICIT NONE
INTEGER(KIND=C_INTPTR_T),INTENT(IN) :: PLAN_LEGACY
TYPE(C_PTR) :: PLAN
PLAN=TRANSFER(PLAN_LEGACY,PLAN)
END FUNCTION TPM_FFTW_CPTR_FROM_LEGACY_PLAN


FUNCTION TPM_FFTW_LEGACY_PLAN_FROM_CPTR(PLAN) RESULT(PLAN_LEGACY)
IMPLICIT NONE
TYPE(C_PTR),INTENT(IN) :: PLAN
INTEGER(KIND=C_INTPTR_T) :: PLAN_LEGACY
PLAN_LEGACY=TRANSFER(PLAN,PLAN_LEGACY)
END FUNCTION TPM_FFTW_LEGACY_PLAN_FROM_CPTR


FUNCTION TPM_FFTW_PLAN_MANY_DFT_C2R_LEGACY(RANK,N,HOWMANY,IN,INEMBED,ISTRIDE,IDIST,OUT,ONEMBED,OSTRIDE,ODIST,FLAGS) RESULT(PLAN)
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN) :: RANK, HOWMANY
INTEGER(KIND=JPIM),INTENT(IN) :: ISTRIDE, IDIST, OSTRIDE, ODIST, FLAGS
INTEGER(KIND=JPIM),INTENT(IN) :: N(*), INEMBED(*), ONEMBED(*)
COMPLEX(KIND=C_DOUBLE_COMPLEX),INTENT(INOUT) :: IN(*)
REAL(KIND=C_DOUBLE),INTENT(INOUT) :: OUT(*)
TYPE(C_PTR) :: PLAN
INTEGER(KIND=C_INTPTR_T) :: PLAN_LEGACY

CALL DFFTW_PLAN_MANY_DFT_C2R(PLAN_LEGACY,RANK,N,HOWMANY,IN,INEMBED,ISTRIDE,IDIST,&
  & OUT,ONEMBED,OSTRIDE,ODIST,FLAGS)
PLAN=TPM_FFTW_CPTR_FROM_LEGACY_PLAN(PLAN_LEGACY)
END FUNCTION TPM_FFTW_PLAN_MANY_DFT_C2R_LEGACY


FUNCTION TPM_FFTWF_PLAN_MANY_DFT_C2R_LEGACY(RANK,N,HOWMANY,IN,INEMBED,ISTRIDE,IDIST,OUT,ONEMBED,OSTRIDE,ODIST,FLAGS) RESULT(PLAN)
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN) :: RANK, HOWMANY
INTEGER(KIND=JPIM),INTENT(IN) :: ISTRIDE, IDIST, OSTRIDE, ODIST, FLAGS
INTEGER(KIND=JPIM),INTENT(IN) :: N(*), INEMBED(*), ONEMBED(*)
COMPLEX(KIND=C_FLOAT_COMPLEX),INTENT(INOUT) :: IN(*)
REAL(KIND=C_FLOAT),INTENT(INOUT) :: OUT(*)
TYPE(C_PTR) :: PLAN
INTEGER(KIND=C_INTPTR_T) :: PLAN_LEGACY

CALL SFFTW_PLAN_MANY_DFT_C2R(PLAN_LEGACY,RANK,N,HOWMANY,IN,INEMBED,ISTRIDE,IDIST,&
  & OUT,ONEMBED,OSTRIDE,ODIST,FLAGS)
PLAN=TPM_FFTW_CPTR_FROM_LEGACY_PLAN(PLAN_LEGACY)
END FUNCTION TPM_FFTWF_PLAN_MANY_DFT_C2R_LEGACY


FUNCTION TPM_FFTW_PLAN_MANY_DFT_R2C_LEGACY(RANK,N,HOWMANY,IN,INEMBED,ISTRIDE,IDIST,OUT,ONEMBED,OSTRIDE,ODIST,FLAGS) RESULT(PLAN)
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN) :: RANK, HOWMANY
INTEGER(KIND=JPIM),INTENT(IN) :: ISTRIDE, IDIST, OSTRIDE, ODIST, FLAGS
INTEGER(KIND=JPIM),INTENT(IN) :: N(*), INEMBED(*), ONEMBED(*)
REAL(KIND=C_DOUBLE),INTENT(INOUT) :: IN(*)
COMPLEX(KIND=C_DOUBLE_COMPLEX),INTENT(INOUT) :: OUT(*)
TYPE(C_PTR) :: PLAN
INTEGER(KIND=C_INTPTR_T) :: PLAN_LEGACY

CALL DFFTW_PLAN_MANY_DFT_R2C(PLAN_LEGACY,RANK,N,HOWMANY,IN,INEMBED,ISTRIDE,IDIST,&
  & OUT,ONEMBED,OSTRIDE,ODIST,FLAGS)
PLAN=TPM_FFTW_CPTR_FROM_LEGACY_PLAN(PLAN_LEGACY)
END FUNCTION TPM_FFTW_PLAN_MANY_DFT_R2C_LEGACY


FUNCTION TPM_FFTWF_PLAN_MANY_DFT_R2C_LEGACY(RANK,N,HOWMANY,IN,INEMBED,ISTRIDE,IDIST,OUT,ONEMBED,OSTRIDE,ODIST,FLAGS) RESULT(PLAN)
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN) :: RANK, HOWMANY
INTEGER(KIND=JPIM),INTENT(IN) :: ISTRIDE, IDIST, OSTRIDE, ODIST, FLAGS
INTEGER(KIND=JPIM),INTENT(IN) :: N(*), INEMBED(*), ONEMBED(*)
REAL(KIND=C_FLOAT),INTENT(INOUT) :: IN(*)
COMPLEX(KIND=C_FLOAT_COMPLEX),INTENT(INOUT) :: OUT(*)
TYPE(C_PTR) :: PLAN
INTEGER(KIND=C_INTPTR_T) :: PLAN_LEGACY

CALL SFFTW_PLAN_MANY_DFT_R2C(PLAN_LEGACY,RANK,N,HOWMANY,IN,INEMBED,ISTRIDE,IDIST,&
  & OUT,ONEMBED,OSTRIDE,ODIST,FLAGS)
PLAN=TPM_FFTW_CPTR_FROM_LEGACY_PLAN(PLAN_LEGACY)
END FUNCTION TPM_FFTWF_PLAN_MANY_DFT_R2C_LEGACY


SUBROUTINE TPM_FFTW_EXECUTE_DFT_C2R_LEGACY(PLAN,IN,OUT)
IMPLICIT NONE
TYPE(C_PTR),INTENT(IN) :: PLAN
COMPLEX(KIND=C_DOUBLE_COMPLEX),INTENT(INOUT) :: IN(*)
REAL(KIND=C_DOUBLE),INTENT(INOUT) :: OUT(*)
INTEGER(KIND=C_INTPTR_T) :: PLAN_LEGACY

PLAN_LEGACY=TPM_FFTW_LEGACY_PLAN_FROM_CPTR(PLAN)
CALL DFFTW_EXECUTE_DFT_C2R(PLAN_LEGACY,IN,OUT)
END SUBROUTINE TPM_FFTW_EXECUTE_DFT_C2R_LEGACY


SUBROUTINE TPM_FFTWF_EXECUTE_DFT_C2R_LEGACY(PLAN,IN,OUT)
IMPLICIT NONE
TYPE(C_PTR),INTENT(IN) :: PLAN
COMPLEX(KIND=C_FLOAT_COMPLEX),INTENT(INOUT) :: IN(*)
REAL(KIND=C_FLOAT),INTENT(INOUT) :: OUT(*)
INTEGER(KIND=C_INTPTR_T) :: PLAN_LEGACY

PLAN_LEGACY=TPM_FFTW_LEGACY_PLAN_FROM_CPTR(PLAN)
CALL SFFTW_EXECUTE_DFT_C2R(PLAN_LEGACY,IN,OUT)
END SUBROUTINE TPM_FFTWF_EXECUTE_DFT_C2R_LEGACY


SUBROUTINE TPM_FFTW_EXECUTE_DFT_R2C_LEGACY(PLAN,IN,OUT)
IMPLICIT NONE
TYPE(C_PTR),INTENT(IN) :: PLAN
REAL(KIND=C_DOUBLE),INTENT(INOUT) :: IN(*)
COMPLEX(KIND=C_DOUBLE_COMPLEX),INTENT(INOUT) :: OUT(*)
INTEGER(KIND=C_INTPTR_T) :: PLAN_LEGACY

PLAN_LEGACY=TPM_FFTW_LEGACY_PLAN_FROM_CPTR(PLAN)
CALL DFFTW_EXECUTE_DFT_R2C(PLAN_LEGACY,IN,OUT)
END SUBROUTINE TPM_FFTW_EXECUTE_DFT_R2C_LEGACY


SUBROUTINE TPM_FFTWF_EXECUTE_DFT_R2C_LEGACY(PLAN,IN,OUT)
IMPLICIT NONE
TYPE(C_PTR),INTENT(IN) :: PLAN
REAL(KIND=C_FLOAT),INTENT(INOUT) :: IN(*)
COMPLEX(KIND=C_FLOAT_COMPLEX),INTENT(INOUT) :: OUT(*)
INTEGER(KIND=C_INTPTR_T) :: PLAN_LEGACY

PLAN_LEGACY=TPM_FFTW_LEGACY_PLAN_FROM_CPTR(PLAN)
CALL SFFTW_EXECUTE_DFT_R2C(PLAN_LEGACY,IN,OUT)
END SUBROUTINE TPM_FFTWF_EXECUTE_DFT_R2C_LEGACY
#endif

END MODULE TPM_FFTW
