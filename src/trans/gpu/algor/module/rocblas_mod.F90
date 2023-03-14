MODULE ROCBLAS_MOD
!
! Define the interfaces to the NVIDIA C code
!
interface roc_gemm
!
! void rocblasSgemm (char transa, char transb, int m, int n,
! int k, float alpha, const float *A, int lda,
! const float *B, int ldb, float beta, float *C, int ldc)
!
subroutine roc_sgemm(cta, ctb, m, n, k,&
alpha, A, lda, B, ldb, beta, c, ldc) bind(C,name='rocblasSgemm')
use iso_c_binding
character(1,c_char),value :: cta, ctb
integer(c_int),value :: m,n,k,lda,ldb,ldc
real(c_float),value :: alpha,beta
real(c_float), dimension(lda,*) :: A
real(c_float), dimension(ldb,*) :: B
real(c_float), dimension(ldc,*) :: C
end subroutine roc_sgemm

!
! void rocblasDgemm (char transa, char transb, int m, int n,
! int k, double alpha, const double *A, int lda,
! const double *B, int ldb, double beta, double *C, int ldc)
!
subroutine roc_dgemm(cta, ctb, m, n, k,&
alpha, A, lda, B, ldb, beta, c, ldc) bind(C,name='rocblasDgemm')
use iso_c_binding
character(1,c_char),value :: cta, ctb
integer(c_int),value :: m,n,k,lda,ldb,ldc
real(c_double),value :: alpha,beta
real(c_double), dimension(lda,*) :: A
real(c_double), dimension(ldb,*) :: B
real(c_double), dimension(ldc,*) :: C
end subroutine roc_dgemm
end interface


INTERFACE
    SUBROUTINE ROC_DGEMM_BATCHED(&
        & CTA, CTB,               &
        & M, N, K,                &
        & ALPHA,                  &
        & A, LDA, TDA,            &
        & B, LDB, TDB,            &
        & BETA,                   &
        & C, LDC, TDC,            &
        & BATCHCOUNT              &
    &) BIND(C, NAME='rocblasDgemmBatched_wrapper')
        USE ISO_C_BINDING
        CHARACTER(1,C_CHAR), VALUE            :: CTA, CTB
        INTEGER(C_INT),      VALUE            :: M, N, K, LDA, LDB, LDC, TDA, TDB, TDC, BATCHCOUNT
        REAL(C_DOUBLE),      VALUE            :: ALPHA,BETA
        REAL(C_DOUBLE),      DIMENSION(LDA,*) :: A
        REAL(C_DOUBLE),      DIMENSION(LDB,*) :: B
        REAL(C_DOUBLE),      DIMENSION(LDC,*) :: C
    END SUBROUTINE ROC_DGEMM_BATCHED

    SUBROUTINE ROC_DGEMM_STRIDED_BATCHED(&
        & CTA, CTB,               &
        & M, N, K,                &
        & ALPHA,                  &
        & A, LDA, TDA,            &
        & B, LDB, TDB,            &
        & BETA,                   &
        & C, LDC, TDC,            &
        & BATCHCOUNT              &
    &) BIND(C, NAME='rocblasDgemmStridedBatched_wrapper')
        USE ISO_C_BINDING
        CHARACTER(1,C_CHAR),  VALUE            :: CTA, CTB
        INTEGER(C_INT),       VALUE            :: M, N, K, LDA, LDB, LDC, BATCHCOUNT
        INTEGER(C_LONG_LONG), VALUE            :: TDA,TDB,TDC
        REAL(C_DOUBLE),        VALUE            :: ALPHA, BETA
        REAL(C_DOUBLE),        DIMENSION(LDA,*) :: A
        REAL(C_DOUBLE),        DIMENSION(LDB,*) :: B
        REAL(C_DOUBLE),        DIMENSION(LDC,*) :: C
    END SUBROUTINE ROC_DGEMM_STRIDED_BATCHED

    subroutine roc_dgemm_batched_finalize() bind(C,name='rocblasDgemmBatched_finalize')
    end subroutine roc_dgemm_batched_finalize

END INTERFACE 

INTERFACE

    SUBROUTINE ROC_SGEMM_BATCHED(&
        & CTA, CTB,               &
        & M, N, K,                &
        & ALPHA,                  &
        & A, LDA, TDA,            &
        & B, LDB, TDB,            &
        & BETA,                   &
        & C, LDC, TDC,            &
        & BATCHCOUNT              &
    &) BIND(C, NAME='rocblasSgemmBatched_wrapper')
        USE ISO_C_BINDING
        CHARACTER(1,C_CHAR), VALUE            :: CTA, CTB
        INTEGER(C_INT),      VALUE            :: M, N, K, LDA, LDB, LDC, TDA, TDB, TDC, BATCHCOUNT
        REAL(C_FLOAT),       VALUE            :: ALPHA, BETA
        REAL(C_FLOAT),       DIMENSION(LDA,*) :: A
        REAL(C_FLOAT),       DIMENSION(LDB,*) :: B
        REAL(C_FLOAT),       DIMENSION(LDC,*) :: C
    END SUBROUTINE ROC_SGEMM_BATCHED
!!END INTERFACE

!!INTERFACE
    SUBROUTINE ROC_SGEMM_STRIDED_BATCHED(&
        & CTA, CTB,               &
        & M, N, K,                &
        & ALPHA,                  &
        & A, LDA, TDA,            &
        & B, LDB, TDB,            &
        & BETA,                   &
        & C, LDC, TDC,            &
        & BATCHCOUNT              &
    &) BIND(C, NAME='rocblasSgemmStridedBatched_wrapper')
        USE ISO_C_BINDING
        CHARACTER(1,C_CHAR),  VALUE            :: CTA, CTB
        INTEGER(C_INT),       VALUE            :: M, N, K, LDA, LDB, LDC, BATCHCOUNT
        INTEGER(C_LONG_LONG), VALUE            :: TDA,TDB,TDC
        REAL(C_FLOAT),        VALUE            :: ALPHA, BETA
        REAL(C_FLOAT),        DIMENSION(LDA,*) :: A
        REAL(C_FLOAT),        DIMENSION(LDB,*) :: B
        REAL(C_FLOAT),        DIMENSION(LDC,*) :: C
    END SUBROUTINE ROC_SGEMM_STRIDED_BATCHED

    subroutine roc_sgemm_batched_finalize() bind(C,name='rocblasSgemmBatched_finalize')
    end subroutine roc_sgemm_batched_finalize


END INTERFACE

INTERFACE
    SUBROUTINE ROC_STCGEMM_BATCHED(&
        & CTA, CTB,                &
        & M, N, K,                 &
        & ALPHA,                   &
        & A, LDA, TDA,             &
        & B, LDB, TDB,             &
        & BETA,                    &
        & C, LDC, TDC,             &
        & BATCHCOUNT               &
    &) BIND(C, NAME='rocblasSTCgemmBatched_wrapper')
        USE ISO_C_BINDING
        CHARACTER(1,C_CHAR), VALUE            :: CTA, CTB
        INTEGER(C_INT),      VALUE            :: M, N, K, LDA, LDB, LDC, TDA, TDB, TDC, BATCHCOUNT
        REAL(C_FLOAT),       VALUE            :: ALPHA, BETA
        REAL(C_FLOAT),       DIMENSION(LDA,*) :: A
        REAL(C_FLOAT),       DIMENSION(LDB,*) :: B
        REAL(C_FLOAT),       DIMENSION(LDC,*) :: C
    END SUBROUTINE ROC_STCGEMM_BATCHED
END INTERFACE




END MODULE ROCBLAS_MOD
