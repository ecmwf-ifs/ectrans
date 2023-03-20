MODULE HIPBLAS_MOD
!
! Define the interfaces to the NVIDIA C code
!
interface hip_gemm
!
! void hipblasSgemm (char transa, char transb, int m, int n,
! int k, float alpha, const float *A, int lda,
! const float *B, int ldb, float beta, float *C, int ldc)
!
subroutine hip_sgemm(cta, ctb, m, n, k,&
alpha, A, lda, B, ldb, beta, c, ldc) bind(C,name='hipblasSgemm')
use iso_c_binding
character(1,c_char),value :: cta, ctb
integer(c_int),value :: m,n,k,lda,ldb,ldc
real(c_float),value :: alpha,beta
real(c_float), dimension(lda,*) :: A
real(c_float), dimension(ldb,*) :: B
real(c_float), dimension(ldc,*) :: C
end subroutine hip_sgemm

!
! void hipblasDgemm (char transa, char transb, int m, int n,
! int k, double alpha, const double *A, int lda,
! const double *B, int ldb, double beta, double *C, int ldc)
!
subroutine hip_dgemm(cta, ctb, m, n, k,&
alpha, A, lda, B, ldb, beta, c, ldc) bind(C,name='hipblasDgemm')
use iso_c_binding
character(1,c_char),value :: cta, ctb
integer(c_int),value :: m,n,k,lda,ldb,ldc
real(c_double),value :: alpha,beta
real(c_double), dimension(lda,*) :: A
real(c_double), dimension(ldb,*) :: B
real(c_double), dimension(ldc,*) :: C
end subroutine hip_dgemm
end interface


INTERFACE
    SUBROUTINE HIP_DGEMM_BATCHED(&
        & CTA, CTB,               &
        & M, N, K,                &
        & ALPHA,                  &
        & A, LDA, TDA,            &
        & B, LDB, TDB,            &
        & BETA,                   &
        & C, LDC, TDC,            &
        & BATCHCOUNT              &
    &) BIND(C, NAME='hipblasDgemmBatched_wrapper')
        USE ISO_C_BINDING
        CHARACTER(1,C_CHAR), VALUE            :: CTA, CTB
        INTEGER(C_INT),      VALUE            :: M, N, K, LDA, LDB, LDC, TDA, TDB, TDC, BATCHCOUNT
        REAL(C_DOUBLE),      VALUE            :: ALPHA,BETA
        REAL(C_DOUBLE),      DIMENSION(LDA,*) :: A
        REAL(C_DOUBLE),      DIMENSION(LDB,*) :: B
        REAL(C_DOUBLE),      DIMENSION(LDC,*) :: C
    END SUBROUTINE HIP_DGEMM_BATCHED

    SUBROUTINE HIP_DGEMM_STRIDED_BATCHED(&
        & CTA, CTB,               &
        & M, N, K,                &
        & ALPHA,                  &
        & A, LDA, TDA,            &
        & B, LDB, TDB,            &
        & BETA,                   &
        & C, LDC, TDC,            &
        & BATCHCOUNT              &
    &) BIND(C, NAME='hipblasDgemmStridedBatched_wrapper')
        USE ISO_C_BINDING
        CHARACTER(1,C_CHAR),  VALUE            :: CTA, CTB
        INTEGER(C_INT),       VALUE            :: M, N, K, LDA, LDB, LDC, BATCHCOUNT
        INTEGER(C_LONG_LONG), VALUE            :: TDA,TDB,TDC
        REAL(C_DOUBLE),        VALUE            :: ALPHA, BETA
        REAL(C_DOUBLE),        DIMENSION(LDA,*) :: A
        REAL(C_DOUBLE),        DIMENSION(LDB,*) :: B
        REAL(C_DOUBLE),        DIMENSION(LDC,*) :: C
    END SUBROUTINE HIP_DGEMM_STRIDED_BATCHED

    subroutine hip_dgemm_batched_finalize() bind(C,name='hipblasDgemmBatched_finalize')
    end subroutine hip_dgemm_batched_finalize

END INTERFACE 

INTERFACE

    SUBROUTINE HIP_SGEMM_BATCHED(&
        & CTA, CTB,               &
        & M, N, K,                &
        & ALPHA,                  &
        & A, LDA, TDA,            &
        & B, LDB, TDB,            &
        & BETA,                   &
        & C, LDC, TDC,            &
        & BATCHCOUNT              &
    &) BIND(C, NAME='hipblasSgemmBatched_wrapper')
        USE ISO_C_BINDING
        CHARACTER(1,C_CHAR), VALUE            :: CTA, CTB
        INTEGER(C_INT),      VALUE            :: M, N, K, LDA, LDB, LDC, TDA, TDB, TDC, BATCHCOUNT
        REAL(C_FLOAT),       VALUE            :: ALPHA, BETA
        REAL(C_FLOAT),       DIMENSION(LDA,*) :: A
        REAL(C_FLOAT),       DIMENSION(LDB,*) :: B
        REAL(C_FLOAT),       DIMENSION(LDC,*) :: C
    END SUBROUTINE HIP_SGEMM_BATCHED
!!END INTERFACE

!!INTERFACE
    SUBROUTINE HIP_SGEMM_STRIDED_BATCHED(&
        & CTA, CTB,               &
        & M, N, K,                &
        & ALPHA,                  &
        & A, LDA, TDA,            &
        & B, LDB, TDB,            &
        & BETA,                   &
        & C, LDC, TDC,            &
        & BATCHCOUNT              &
    &) BIND(C, NAME='hipblasSgemmStridedBatched_wrapper')
        USE ISO_C_BINDING
        CHARACTER(1,C_CHAR),  VALUE            :: CTA, CTB
        INTEGER(C_INT),       VALUE            :: M, N, K, LDA, LDB, LDC, BATCHCOUNT
        INTEGER(C_LONG_LONG), VALUE            :: TDA,TDB,TDC
        REAL(C_FLOAT),        VALUE            :: ALPHA, BETA
        REAL(C_FLOAT),        DIMENSION(LDA,*) :: A
        REAL(C_FLOAT),        DIMENSION(LDB,*) :: B
        REAL(C_FLOAT),        DIMENSION(LDC,*) :: C
    END SUBROUTINE HIP_SGEMM_STRIDED_BATCHED

    subroutine hip_sgemm_batched_finalize() bind(C,name='hipblasSgemmBatched_finalize')
    end subroutine hip_sgemm_batched_finalize


END INTERFACE

INTERFACE
    SUBROUTINE HIP_STCGEMM_BATCHED(&
        & CTA, CTB,                &
        & M, N, K,                 &
        & ALPHA,                   &
        & A, LDA, TDA,             &
        & B, LDB, TDB,             &
        & BETA,                    &
        & C, LDC, TDC,             &
        & BATCHCOUNT               &
    &) BIND(C, NAME='hipblasSTCgemmBatched_wrapper')
        USE ISO_C_BINDING
        CHARACTER(1,C_CHAR), VALUE            :: CTA, CTB
        INTEGER(C_INT),      VALUE            :: M, N, K, LDA, LDB, LDC, TDA, TDB, TDC, BATCHCOUNT
        REAL(C_FLOAT),       VALUE            :: ALPHA, BETA
        REAL(C_FLOAT),       DIMENSION(LDA,*) :: A
        REAL(C_FLOAT),       DIMENSION(LDB,*) :: B
        REAL(C_FLOAT),       DIMENSION(LDC,*) :: C
    END SUBROUTINE HIP_STCGEMM_BATCHED
END INTERFACE




END MODULE HIPBLAS_MOD
