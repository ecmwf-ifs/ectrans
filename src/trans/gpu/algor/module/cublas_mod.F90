MODULE CUBLAS_MOD
!
! Define the interfaces to the NVIDIA C code
!
interface cuda_gemm
!
! void cublasSgemm (char transa, char transb, int m, int n,
! int k, float alpha, const float *A, int lda,
! const float *B, int ldb, float beta, float *C, int ldc)
!
subroutine cuda_sgemm(cta, ctb, m, n, k,&
alpha, A, lda, B, ldb, beta, c, ldc) bind(C,name='cublasSgemm')
use iso_c_binding
character(1,c_char),value :: cta, ctb
integer(c_int),value :: m,n,k,lda,ldb,ldc
real(c_float),value :: alpha,beta
real(c_float), dimension(lda,*) :: A
real(c_float), dimension(ldb,*) :: B
real(c_float), dimension(ldc,*) :: C
end subroutine cuda_sgemm

!
! void cublasDgemm (char transa, char transb, int m, int n,
! int k, double alpha, const double *A, int lda,
! const double *B, int ldb, double beta, double *C, int ldc)
!
subroutine cuda_dgemm(cta, ctb, m, n, k,&
alpha, A, lda, B, ldb, beta, c, ldc) bind(C,name='cublasDgemm')
use iso_c_binding
character(1,c_char),value :: cta, ctb
integer(c_int),value :: m,n,k,lda,ldb,ldc
real(c_double),value :: alpha,beta
real(c_double), dimension(lda,*) :: A
real(c_double), dimension(ldb,*) :: B
real(c_double), dimension(ldc,*) :: C
end subroutine cuda_dgemm
end interface


INTERFACE
    SUBROUTINE CUDA_DGEMM_BATCHED(&
        & CTA, CTB,               &
        & M, N, K,                &
        & ALPHA,                  &
        & A, LDA, TDA,            &
        & B, LDB, TDB,            &
        & BETA,                   &
        & C, LDC, TDC,            &
        & BATCHCOUNT              &
    &) BIND(C, NAME='cublasDgemmBatched_wrapper')
        USE ISO_C_BINDING
        CHARACTER(1,C_CHAR), VALUE            :: CTA, CTB
        INTEGER(C_INT),      VALUE            :: M, N, K, LDA, LDB, LDC, TDA, TDB, TDC, BATCHCOUNT
        REAL(C_DOUBLE),      VALUE            :: ALPHA,BETA
        REAL(C_DOUBLE),      DIMENSION(LDA,*) :: A
        REAL(C_DOUBLE),      DIMENSION(LDB,*) :: B
        REAL(C_DOUBLE),      DIMENSION(LDC,*) :: C
    END SUBROUTINE CUDA_DGEMM_BATCHED

    SUBROUTINE CUDA_DGEMM_STRIDED_BATCHED(&
        & CTA, CTB,               &
        & M, N, K,                &
        & ALPHA,                  &
        & A, LDA, TDA,            &
        & B, LDB, TDB,            &
        & BETA,                   &
        & C, LDC, TDC,            &
        & BATCHCOUNT              &
    &) BIND(C, NAME='cublasDgemmStridedBatched_wrapper')
        USE ISO_C_BINDING
        CHARACTER(1,C_CHAR),  VALUE            :: CTA, CTB
        INTEGER(C_INT),       VALUE            :: M, N, K, LDA, LDB, LDC, BATCHCOUNT
        INTEGER(C_LONG_LONG), VALUE            :: TDA,TDB,TDC
        REAL(C_DOUBLE),        VALUE            :: ALPHA, BETA
        REAL(C_DOUBLE),        DIMENSION(LDA,*) :: A
        REAL(C_DOUBLE),        DIMENSION(LDB,*) :: B
        REAL(C_DOUBLE),        DIMENSION(LDC,*) :: C
    END SUBROUTINE CUDA_DGEMM_STRIDED_BATCHED

    subroutine cuda_dgemm_batched_finalize() bind(C,name='cublasDgemmBatched_finalize')
    end subroutine cuda_dgemm_batched_finalize

END INTERFACE 

INTERFACE

    SUBROUTINE CUDA_SGEMM_BATCHED(&
        & CTA, CTB,               &
        & M, N, K,                &
        & ALPHA,                  &
        & A, LDA, TDA,            &
        & B, LDB, TDB,            &
        & BETA,                   &
        & C, LDC, TDC,            &
        & BATCHCOUNT              &
    &) BIND(C, NAME='cublasSgemmBatched_wrapper')
        USE ISO_C_BINDING
        CHARACTER(1,C_CHAR), VALUE            :: CTA, CTB
        INTEGER(C_INT),      VALUE            :: M, N, K, LDA, LDB, LDC, TDA, TDB, TDC, BATCHCOUNT
        REAL(C_FLOAT),       VALUE            :: ALPHA, BETA
        REAL(C_FLOAT),       DIMENSION(LDA,*) :: A
        REAL(C_FLOAT),       DIMENSION(LDB,*) :: B
        REAL(C_FLOAT),       DIMENSION(LDC,*) :: C
    END SUBROUTINE CUDA_SGEMM_BATCHED
!!END INTERFACE

!!INTERFACE
    SUBROUTINE CUDA_SGEMM_STRIDED_BATCHED(&
        & CTA, CTB,               &
        & M, N, K,                &
        & ALPHA,                  &
        & A, LDA, TDA,            &
        & B, LDB, TDB,            &
        & BETA,                   &
        & C, LDC, TDC,            &
        & BATCHCOUNT              &
    &) BIND(C, NAME='cublasSgemmStridedBatched_wrapper')
        USE ISO_C_BINDING
        CHARACTER(1,C_CHAR),  VALUE            :: CTA, CTB
        INTEGER(C_INT),       VALUE            :: M, N, K, LDA, LDB, LDC, BATCHCOUNT
        INTEGER(C_LONG_LONG), VALUE            :: TDA,TDB,TDC
        REAL(C_FLOAT),        VALUE            :: ALPHA, BETA
        REAL(C_FLOAT),        DIMENSION(LDA,*) :: A
        REAL(C_FLOAT),        DIMENSION(LDB,*) :: B
        REAL(C_FLOAT),        DIMENSION(LDC,*) :: C
    END SUBROUTINE CUDA_SGEMM_STRIDED_BATCHED

    subroutine cuda_sgemm_batched_finalize() bind(C,name='cublasSgemmBatched_finalize')
    end subroutine cuda_sgemm_batched_finalize


END INTERFACE

INTERFACE
    SUBROUTINE CUDA_STCGEMM_BATCHED(&
        & CTA, CTB,                &
        & M, N, K,                 &
        & ALPHA,                   &
        & A, LDA, TDA,             &
        & B, LDB, TDB,             &
        & BETA,                    &
        & C, LDC, TDC,             &
        & BATCHCOUNT               &
    &) BIND(C, NAME='cublasSTCgemmBatched_wrapper')
        USE ISO_C_BINDING
        CHARACTER(1,C_CHAR), VALUE            :: CTA, CTB
        INTEGER(C_INT),      VALUE            :: M, N, K, LDA, LDB, LDC, TDA, TDB, TDC, BATCHCOUNT
        REAL(C_FLOAT),       VALUE            :: ALPHA, BETA
        REAL(2),             DIMENSION(LDA,*) :: A
        REAL(2),             DIMENSION(LDB,*) :: B
        REAL(C_FLOAT),       DIMENSION(LDC,*) :: C
    END SUBROUTINE CUDA_STCGEMM_BATCHED
END INTERFACE




END MODULE CUBLAS_MOD
