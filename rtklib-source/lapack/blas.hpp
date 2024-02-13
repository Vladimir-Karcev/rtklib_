#ifndef BLAS_HPP_MAGS
#define BLAS_HPP_MAGS

#include "f2c.h"

#undef abs
#undef dabs
#undef min
#undef max
#undef dmin
#undef dmax
#undef bit_test
#undef bit_clear
#undef bit_set

#include "blaswrap.h"
#include "clapack.h"

// ASUM takes the sum of the absolute values.
inline double asum(long n, const double *x, long incx) {
	return dasum_(&n, (double *)x, &incx);
}

// AXPY constant times a vector plus a vector.
inline void axpy(long n, double alpha, const double *x, long incx, double *y,
	long incy) {
	daxpy_(&n, &alpha, (double *)x, &incx, y, &incy);
}

// COPY copies a vector, x, to a vector, y.
inline void copy(long n, const double *x, long incx, double *y, long incy) {
	dcopy_(&n, (double *)x, &incx, y, &incy);
}

// DOT forms the dot product of two vectors.
inline double dot(long n, const double *x, long incx, const double *y, long incy) {
	return ddot_(&n, (double *)x, &incx, (double *)y, &incy);
}

// NRM2 returns the euclidean norm of a vector.
inline double nrm2(long n, const double *x, long incx) {
	return dnrm2_(&n, (double *)x, &incx);
}

// SCAL scales a vector by a constant.
inline void scal(long n, double a, double*x, long incx) {
	dscal_(&n, &a, x, &incx);
}

// SWAP interchanges two vectors.
inline void swap(long n, double* x, long incx, double* y, long incy) {
	dswap_(&n, x, &incx, y, &incy);
}

// GEMV performs one of the matrix-vector operations
// 
//    y : = alpha*A*x + beta*y, or   y : = alpha*A**T*x + beta*y,
// 
// where alpha and beta are scalars, x and y are vectors and A is an m by n
// matrix.
inline void gemv(const char *trans, long m, long n, double alpha,
	const double *a, long lda, const double *x, long incx, double beta,
	double *y, long incy) {
	dgemv_((char *)trans, &m, &n, &alpha, (double *)a, &lda, (double *)x, &incx,
		&beta, y, &incy);
}

// GER performs the rank 1 operation
//
//    A : = alpha*x*y**T + A,
//
// where alpha is a scalar, x is an m element vector, y is an n element vector
// and A is an m by n matrix.
inline void ger(long m, long n, double alpha, const double *x, long incx,
	const double *y, long incy, double *a, long lda) {
	dger_(&m, &n, &alpha, (double *)x, &incx, (double *)y, &incy, a, &lda);
}

// TRMV performs one of the matrix-vector operations
//
//    x : = A*x, or x : = A**T*x,
//
// where x is an n element vector and  A is an n by n unit, or non-unit, upper
// or lower triangular matrix.
inline void trmv(const char *uplo, const char *trans, const char *diag, long n,
	const double *a, long lda, double *x, long incx) {
	dtrmv_((char *)uplo, (char *)trans, (char *)diag, &n, (double *)a, &lda, x, &incx);
}

// DTRSM  solves one of the matrix equations
//
//    op(A) *X = alpha*B, or   X*op(A) = alpha*B,
//
// where alpha is a scalar, X andB are m by n matrices, A is a unit, or
// non-unit, upper or lower triangular matrix  and  op(A)  is one  of
//
//    op(A) = A   or   op(A) = A**T.
//
// The matrix X is overwritten on B.
inline void trsm(const char *side, const char *uplo, const char *trans, const char *diag,
	long m, long n, double alpha, const double *a, long lda, double *b, long ldb) {
	dtrsm_((char *)side, (char *)uplo, (char *)trans, (char *)diag, &m, &n, &alpha, (double *)a, &lda, b, &ldb);
}

// GEMM performs one of the matrix-matrix operations
// 
//    C : = alpha*op(A)*op(B) + beta*C,
// 
// where  op(X) is one of
// 
//    op(X) = X   or   op(X) = X**T,
// 
// alpha and beta are scalars, and A, B and C are matrices, with op(A) an
// m by k matrix, op(B) a k by n matrix and C an m by n matrix.
inline void gemm(const char *transa, const char *transb, long m, long n,
	long k, double alpha, const double *a, long lda, const double *b, long ldb,
	double beta, double *c, long ldc) {
	dgemm_((char *)transa, (char *)transb, &m, &n, &k, &alpha, (double *)a,
		&lda, (double *)b, &ldb, &beta, c, &ldc);
}

// SYMV performs the matrix-vector operation
//
//    y := alpha*A*x + beta*y,
//
//  where alpha and beta are scalars, x and y are n element vectors and A is
//  an n by n symmetric matrix.
inline void symv(const char *uplo, long n, double alpha, const double *a, long lda,
	const double *x, long incx, double beta, double *y, long incy) {
	dsymv_((char *)uplo, &n, &alpha, (double *)a, &lda, (double *)x, &incx,
		&beta, y, &incy);
}

// SYMM  performs one of the matrix-matrix operations
//
//    C := alpha*A*B + beta*C,
//
// or
//
//    C := alpha*B*A + beta*C,
//
// where alpha and beta are scalars, A is a symmetric matrix and B and C are
// m by n matrices.
inline void symm(const char *side, const char *uplo, long m, long n,
	double alpha, const double *a, long lda, const double *b, long ldb,
	double beta, double *c, long ldc) {
	dsymm_((char *)side, (char *)uplo, &m, &n, &alpha, (double *)a,
		&lda, (double *)b, &ldb, &beta, c, &ldc);
}

// SYRK performs one of the symmetric rank k operations
//
//    C : = alpha*A*A**T + beta*C,
//
// or
//
//    C : = alpha*A**T*A + beta*C,
//
// where alpha andbeta are scalars, C is an n by n symmetric matrix and A is an
// n by k matrix in the first case anda k by n matrix in the second case.
inline void syrk(const char *uplo, const char *trans, long n, long k,
	double alpha, const double *a, long lda, double beta, double *c, long ldc) {
	dsyrk_((char *)uplo, (char *)trans, &n, &k, &alpha, (double *)a, &lda, &beta, c, &ldc);
}

#endif
