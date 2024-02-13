#ifndef LAPACK_HPP_MAGS
#define LAPACK_HPP_MAGS

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

#include "clapack.h"

// GETRF computes an LU factorization of a general M-by-N matrix A
// using partial pivoting with row interchanges.
//
// The factorization has the form
//    A = P*L*U
// where P is a permutation matrix, L is lower triangular with unit
// diagonal elements(lower trapezoidal if m > n), andU is upper
// triangular(upper trapezoidal if m < n).
inline long getrf(long m, long n, double *a, long lda, long *ipiv) {
	long info;
	dgetrf_(&m, &n, a, &lda, ipiv, &info);
	return info;
}

// GETRS solves a system of linear equations
//    A*X = B  or  A**T*X = B
// with a general N-by-N matrix A using the LU factorization computed
// by GETRF.
inline long getrs(const char *trans, long n, long nrhs, const double *a, long lda,
	long *ipiv, double *b, long ldb) {
	long info;
	dgetrs_((char *)trans, &n, &nrhs, (double *)a, &lda, ipiv, b, &ldb, &info);
	return info;
}

// GESV computes the solution to a real system of linear equations
//    A*X = B,
// where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
//
// The LU decomposition with partial pivoting androw interchanges is used to
// factor A as
//    A = P*L*U,
// where P is a permutation matrix, L is unit lower triangular, and U is upper
// triangular. The factored form of A is then used to solve the system of
// equations A*X = B.
inline long gesv(long n, long nrhs, const double *a, long lda, long *ipiv,
	double *b, long ldb) {
	long info;
	dgesv_(&n, &nrhs, (double *)a, &lda, ipiv, b, &ldb, &info);
	return info;
}

// TRTRS solves a triangular system of the form
//
//    A *X = B  or  A**T * X = B,
//
// where A is a triangular matrix of order N, and B is an N-by-NRHS matrix.
// A check is made to verify that A is nonsingular.
inline long trtrs(const char *uplo, const char *trans, const char *diag, long n,
	long nrhs, const double *a, long lda, double *b, long ldb) {
	long info;
	dtrtrs_((char *)uplo, (char *)trans, (char *)diag, &n, &nrhs, (double *)a,
		&lda, b, &ldb, &info);
	return info;
}

// POTRF computes the Cholesky factorization of a real symmetric positive
// definite matrix A.
// 
// The factorization has the form
//    A = U**T*U, if UPLO = 'U', or
//    A = L*L**T, if UPLO = 'L',
// where U is an upper triangular matrix and L is lower triangular.
inline long potrf(const char *uplo, long n, double *a, long lda) {
	long info;
	dpotrf_((char *)uplo, &n, a, &lda, &info);
	return info;
}

// POCON estimates the reciprocal of the condition number(in the 1-norm) of
// a real symmetric positive definite matrix using the Cholesky factorization
// A = U**T*U or A = L*L**T computed by POTRF.
// 
// An estimate is obtained for norm(inv(A)), andthe reciprocal of the condition
// number is computed as RCOND = 1/(ANORM*norm(inv(A))).
inline long pocon(const char *uplo, long n, double *a, long lda, double anorm,
	double &rcond, double *work, long *iwork) {
	long info;
	dpocon_((char *)uplo, &n, a, &lda, &anorm, &rcond, work, iwork, &info);
	return info;
}

// POTRS solves a system of linear equations A*X = B with a symmetric positive
// definite matrix A using the Cholesky factorization A = U**T*U or A = L*L**T
// computed by POTRF.
inline long potrs(const char *uplo, long n, long nrhs, const double *a,
	long lda, double *b, long ldb) {
	long info;
	dpotrs_((char *)uplo, &n, &nrhs, (double *)a, &lda, b, &ldb, &info);
	return info;
}

// POTRI computes the inverse of a real symmetric positive definite matrix
// A using the Cholesky factorization A = U**T*U or A = L*L**T computed by
// POTRF.
inline long potri(const char *uplo, long n, double *a, long lda) {
	long info;
	dpotri_((char *)uplo, &n, a, &lda, &info);
	return info;
}

// POSV computes the solution to a real system of linear equations
//    A *X = B,
// where A is an N-by-N symmetric positive definite matrix and X and B are
// N-by-NRHS matrices.
// 
// The Cholesky decomposition is used to factor A as
//    A = U**T* U, if UPLO = 'U', or
//    A = L * L**T, if UPLO = 'L',
// where U is an upper triangular matrix andL is a lower triangular matrix. The
// factored form of A is then used to solve the system of equations A *X = B.
inline long posv(const char *uplo, long n, long nrhs, double *a, long lda,
	double *b, long ldb) {
	long info;
	dposv_((char *)uplo, &n, &nrhs, a, &lda, b, &ldb, &info);
	return info;
}

// POSVX uses the Cholesky factorization A = U**T*U or A = L*L**T to compute
// the solution to a real system of linear equations
//    A *X = B,
// where A is an N-by-N symmetric positive definite matrix and X and B are
// N-by-NRHS matrices.
// 
// Error bounds on the solution anda condition estimate are also provided.
inline long posvx(const char *fact, const char *uplo, long n, long nrhs,
	double *a, long lda, double *af, long ldaf, char *equed, double *s, double *b,
	long ldb, double *x, long ldx, double &rcond, double *ferr, double *berr,
	double *work, long *iwork) {
	long info;
	dposvx_((char *)fact, (char *)uplo, &n, &nrhs, a, &lda, af, &ldaf, equed,
		s, b, &ldb, x, &ldx, &rcond, ferr, berr, work, iwork,
		&info);
	return info;
}

//// SYEVR computes selected eigenvalues and, optionally, eigenvectors of a real
//// symmetric matrix A.Eigenvalues andeigenvectors can be selected by specifying
//// either a range of values or a range of indices for the desired eigenvalues.
//inline long syevr() {
//	long info = 0;
//	return info;
//}
//
//// SYEVD computes all eigenvalues and, optionally, eigenvectors of a real
//// symmetric matrix A.If eigenvectors are desired, it uses a divide and conquer
//// algorithm.
//inline long syevd(const char *jobz, const char *uplo, long n, double *a, long lda,
//	double *w, double *work, long lwork, long *iwork, long liwork) {
//	long info;
//	dsyevd_((char *)jobz, (char *)uplo, &n, a, &lda, w, work, &lwork, iwork, &liwork, &info);
//	return info;
//}

// LACPY copies all or part of a two-dimensional matrix A to another matrix B.
inline void lacpy(const char *uplo, long m, long n, const double *a, long lda,
	double *b, long ldb) {
	dlacpy_((char *)uplo, &m, &n, (double *)a, &lda, b, &ldb);
}

// LANSY returns the value of the one norm, or the Frobenius norm, or the
// infinity norm, or the element of largest absolute value of a real symmetric
// matrix A.
inline double lansy(const char *norm, const char *uplo, long n, const double *a,
	long lda, double *work) {
	return dlansy_((char *)norm, (char *)uplo, &n, (double *)a, &lda, work);
}

// LASET initializes an m-by-n matrix A to BETA on the diagonal and ALPHA
// on the offdiagonals.
inline void laset(const char *uplo, long m, long n, double alpha, double beta,
	double *a, long lda) {
	dlaset_((char *)uplo, &m, &n, &alpha, &beta, a, &lda);
}

#endif
