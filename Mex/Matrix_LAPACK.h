#ifndef Matrix_LAPACK_h
#define Matrix_LAPACK_h

#include <valarray>
#include <complex>
#define MKL_Complex8 std::complex<float>
#define MKL_Complex16 std::complex<double>
#include "a_matrices.h"
#include "mkl_lapacke.h"
#include "Matrix_BLAS.h"
#include "valarray_ComplexRealOperators.h"

using namespace std;
using namespace a_matrices;

inline Matrix<float> Solve(const Matrix<float> & A_mat, const Matrix<float> & B_mat)
{
	int N=A_mat.rows;
	if (N!=A_mat.cols) throw("Matrix A_mat must be square.");
	if (N!=B_mat.rows) throw("Matrix B_mat must have the same number of rows as Matrix A_mat.");

	Matrix<float> A_temp(A_mat), C(B_mat);
	valarray<lapack_int> ipiv(N);
	lapack_int info=LAPACKE_sgesv(LAPACK_ROW_MAJOR, N, B_mat.cols, &A_temp[0][0], A_mat.cols, &ipiv[0], &C[0][0], B_mat.cols);
	if (info != 0) throw("LAPACKE_sgesv failed.");

	return C;
}

inline Matrix<double> Solve(const Matrix<double> & A_mat, const Matrix<double> & B_mat)
{
	int N=A_mat.rows;
	if (N!=A_mat.cols) throw("Matrix A_mat must be square.");
	if (N!=B_mat.rows) throw("Matrix B_mat must have the same number of rows as Matrix A_mat.");

	Matrix<double> A_temp(A_mat), C(B_mat);
	valarray<lapack_int> ipiv(N);
	lapack_int info=LAPACKE_dgesv(LAPACK_ROW_MAJOR, N, B_mat.cols, &A_temp[0][0], A_mat.cols, &ipiv[0], &C[0][0], B_mat.cols);
	if (info != 0) throw("LAPACKE_dgesv failed.");

	return C;
}

inline Matrix<complex<float> > Solve(const Matrix<complex<float> > & A_mat, const Matrix<complex<float> > & B_mat)
{
	int N=A_mat.rows;
	if (N!=A_mat.cols) throw("Matrix A_mat must be square.");
	if (N!=B_mat.rows) throw("Matrix B_mat must have the same number of rows as Matrix A_mat.");

	Matrix<complex<float> > A_temp(A_mat), C(B_mat);
	valarray<lapack_int> ipiv(N);
	lapack_int info=LAPACKE_cgesv(LAPACK_ROW_MAJOR, N, B_mat.cols, &A_temp[0][0], A_mat.cols, &ipiv[0], &C[0][0], B_mat.cols);
	if (info != 0) throw("LAPACKE_cgesv failed.");

	return C;
}

inline Matrix<complex<double> > Solve(const Matrix<complex<double> > & A_mat, const Matrix<complex<double> > & B_mat)
{
	int N=A_mat.rows;
	if (N!=A_mat.cols) throw("Matrix A_mat must be square.");
	if (N!=B_mat.rows) throw("Matrix B_mat must have the same number of rows as Matrix A_mat.");

	Matrix<complex<double> > A_temp(A_mat), C(B_mat);
	valarray<lapack_int> ipiv(N);
	lapack_int info=LAPACKE_zgesv(LAPACK_ROW_MAJOR, N, B_mat.cols, &A_temp[0][0], A_mat.cols, &ipiv[0], &C[0][0], B_mat.cols);
	if (info != 0) throw("LAPACKE_zgesv failed.");

	return C;
}

inline Matrix<float> SymPosDefSolve(const Matrix<float> & A_mat, const Matrix<float> & B_mat)
{
	int N=A_mat.rows;
	if (N!=A_mat.cols) throw("Matrix A_mat must be square.");
	if (N!=B_mat.rows) throw("Matrix B_mat must have the same number of rows as Matrix A_mat.");

	Matrix<float> A_temp(A_mat), C(B_mat);
	lapack_int info=LAPACKE_sposv(LAPACK_ROW_MAJOR, 'U', N, B_mat.cols, &A_temp[0][0], A_mat.cols, &C[0][0], B_mat.cols);
	if (info > 0)
		throw("Error using LAPACKE_sposv, matrix A_mat is not positive definite.");
	else if (info != 0)
		throw("LAPACKE_sposv failed.");

	return C;
}

inline Matrix<double> SymPosDefSolve(const Matrix<double> & A_mat, const Matrix<double> & B_mat)
{
	int N=A_mat.rows;
	if (N!=A_mat.cols) throw("Matrix A_mat must be square.");
	if (N!=B_mat.rows) throw("Matrix B_mat must have the same number of rows as Matrix A_mat.");

	Matrix<double> A_temp(A_mat), C(B_mat);
	lapack_int info=LAPACKE_dposv(LAPACK_ROW_MAJOR, 'U', N, B_mat.cols, &A_temp[0][0], A_mat.cols, &C[0][0], B_mat.cols);
	if (info > 0)
		throw("Error using LAPACKE_dposv, matrix A_mat is not positive definite.");
	else if (info != 0)
		throw("LAPACKE_dposv failed.");

	return C;
}

inline Matrix<complex<float> > HermPosDefSolve(const Matrix<complex<float> > & A_mat, const Matrix<complex<float> > & B_mat)
{
	int N=A_mat.rows;
	if (N!=A_mat.cols) throw("Matrix A_mat must be square.");
	if (N!=B_mat.rows) throw("Matrix B_mat must have the same number of rows as Matrix A_mat.");

	Matrix<complex<float> > A_temp(A_mat), C(B_mat);
	lapack_int info=LAPACKE_cposv(LAPACK_ROW_MAJOR, 'U', N, B_mat.cols, &A_temp[0][0], A_mat.cols, &C[0][0], B_mat.cols);
	if (info > 0)
		throw("Error using LAPACKE_cposv, matrix A_mat is not positive definite.");
	else if (info != 0)
		throw("LAPACKE_cposv failed.");

	return C;
}

inline Matrix<complex<double> > HermPosDefSolve(const Matrix<complex<double> > & A_mat, const Matrix<complex<double> > & B_mat)
{
	int N=A_mat.rows;
	if (N!=A_mat.cols) throw("Matrix A_mat must be square.");
	if (N!=B_mat.rows) throw("Matrix B_mat must have the same number of rows as Matrix A_mat.");

	Matrix<complex<double> > A_temp(A_mat), C(B_mat);
	lapack_int info=LAPACKE_zposv(LAPACK_ROW_MAJOR, 'U', N, B_mat.cols, &A_temp[0][0], A_mat.cols, &C[0][0], B_mat.cols);
	if (info > 0)
		throw("Error using LAPACKE_zposv, matrix A_mat is not positive definite.");
	else if (info != 0)
		throw("LAPACKE_zposv failed.");

	return C;
}

inline Matrix<float> SymSolve(const Matrix<float> & A_mat, const Matrix<float> & B_mat)
{
	int N=A_mat.rows;
	if (N!=A_mat.cols) throw("Matrix A_mat must be square.");
	if (N!=B_mat.rows) throw("Matrix B_mat must have the same number of rows as Matrix A_mat.");

	Matrix<float> A_temp(A_mat), C(B_mat);
	valarray<lapack_int> ipiv(N);
	lapack_int info=LAPACKE_ssysv(LAPACK_ROW_MAJOR, 'U', N, B_mat.cols, &A_temp[0][0], A_mat.cols, &ipiv[0], &C[0][0], B_mat.cols);
	if (info != 0) throw("LAPACKE_ssysv failed.");

	return C;
}

inline Matrix<double> SymSolve(const Matrix<double> & A_mat, const Matrix<double> & B_mat)
{
	int N=A_mat.rows;
	if (N!=A_mat.cols) throw("Matrix A_mat must be square.");
	if (N!=B_mat.rows) throw("Matrix B_mat must have the same number of rows as Matrix A_mat.");

	Matrix<double> A_temp(A_mat), C(B_mat);
	valarray<lapack_int> ipiv(N);
	lapack_int info=LAPACKE_dsysv(LAPACK_ROW_MAJOR, 'U', N, B_mat.cols, &A_temp[0][0], A_mat.cols, &ipiv[0], &C[0][0], B_mat.cols);
	if (info != 0) throw("LAPACKE_dsysv failed.");

	return C;
}

inline Matrix<complex<float> > SymSolve(const Matrix<complex<float> > & A_mat, const Matrix<complex<float> > & B_mat)
{
	int N=A_mat.rows;
	if (N!=A_mat.cols) throw("Matrix A_mat must be square.");
	if (N!=B_mat.rows) throw("Matrix B_mat must have the same number of rows as Matrix A_mat.");

	Matrix<complex<float> > A_temp(A_mat), C(B_mat);
	valarray<lapack_int> ipiv(N);
	lapack_int info=LAPACKE_csysv(LAPACK_ROW_MAJOR, 'U', N, B_mat.cols, &A_temp[0][0], A_mat.cols, &ipiv[0], &C[0][0], B_mat.cols);
	if (info != 0) throw("LAPACKE_csysv failed.");

	return C;
}

inline Matrix<complex<double> > SymSolve(const Matrix<complex<double> > & A_mat, const Matrix<complex<double> > & B_mat)
{
	int N=A_mat.rows;
	if (N!=A_mat.cols) throw("Matrix A_mat must be square.");
	if (N!=B_mat.rows) throw("Matrix B_mat must have the same number of rows as Matrix A_mat.");

	Matrix<complex<double> > A_temp(A_mat), C(B_mat);
	valarray<lapack_int> ipiv(N);
	lapack_int info=LAPACKE_zsysv(LAPACK_ROW_MAJOR, 'U', N, B_mat.cols, &A_temp[0][0], A_mat.cols, &ipiv[0], &C[0][0], B_mat.cols);
	if (info != 0) throw("LAPACKE_zsysv failed.");

	return C;
}

inline Matrix<complex<float> > HermSolve(const Matrix<complex<float> > & A_mat, const Matrix<complex<float> > & B_mat)
{
	int N=A_mat.rows;
	if (N!=A_mat.cols) throw("Matrix A_mat must be square.");
	if (N!=B_mat.rows) throw("Matrix B_mat must have the same number of rows as Matrix A_mat.");

	Matrix<complex<float> > A_temp(A_mat), C(B_mat);
	valarray<lapack_int> ipiv(N);
	lapack_int info=LAPACKE_chesv(LAPACK_ROW_MAJOR, 'U', N, B_mat.cols, &A_temp[0][0], A_mat.cols, &ipiv[0], &C[0][0], B_mat.cols);
	if (info != 0) throw("LAPACKE_chesv failed.");

	return C;
}

inline Matrix<complex<double> > HermSolve(const Matrix<complex<double> > & A_mat, const Matrix<complex<double> > & B_mat)
{
	int N=A_mat.rows;
	if (N!=A_mat.cols) throw("Matrix A_mat must be square.");
	if (N!=B_mat.rows) throw("Matrix B_mat must have the same number of rows as Matrix A_mat.");

	Matrix<complex<double> > A_temp(A_mat), C(B_mat);
	valarray<lapack_int> ipiv(N);
	lapack_int info=LAPACKE_zhesv(LAPACK_ROW_MAJOR, 'U', N, B_mat.cols, &A_temp[0][0], A_mat.cols, &ipiv[0], &C[0][0], B_mat.cols);
	if (info != 0) throw("LAPACKE_zhesv failed.");

	return C;
}

inline Matrix<float> Inv(const Matrix<float> & A_mat)
{
	int N=A_mat.rows;
	if (N!=A_mat.cols) throw("Matrix A_mat must be square.");

	Matrix<float> A_inv(A_mat);
	valarray<lapack_int> ipiv(N);
	lapack_int info=LAPACKE_sgetrf(LAPACK_ROW_MAJOR, N, N, &A_inv[0][0], A_inv.cols, &ipiv[0]);
	if (info != 0) throw("LAPACKE_sgetrf failed.");
	info=LAPACKE_sgetri(LAPACK_ROW_MAJOR, N, &A_inv[0][0], A_inv.cols, &ipiv[0]);
	if (info != 0) throw("LAPACKE_sgetri failed.");

	return A_inv;
}

inline Matrix<double> Inv(const Matrix<double> & A_mat)
{
	int N=A_mat.rows;
	if (N!=A_mat.cols) throw("Matrix A_mat must be square.");

	Matrix<double> A_inv(A_mat);
	valarray<lapack_int> ipiv(N);
	lapack_int info=LAPACKE_dgetrf(LAPACK_ROW_MAJOR, N, N, &A_inv[0][0], A_inv.cols, &ipiv[0]);
	if (info != 0) throw("LAPACKE_dgetrf failed.");
	info=LAPACKE_dgetri(LAPACK_ROW_MAJOR, N, &A_inv[0][0], A_inv.cols, &ipiv[0]);
	if (info != 0) throw("LAPACKE_dgetri failed.");

	return A_inv;
}

inline Matrix<complex<float> > Inv(const Matrix<complex<float> > & A_mat)
{
	int N=A_mat.rows;
	if (N!=A_mat.cols) throw("Matrix A_mat must be square.");

	Matrix<complex<float> > A_inv(A_mat);
	valarray<lapack_int> ipiv(N);
	lapack_int info=LAPACKE_cgetrf(LAPACK_ROW_MAJOR, N, N, &A_inv[0][0], A_inv.cols, &ipiv[0]);
	if (info != 0) throw("LAPACKE_cgetrf failed.");
	info=LAPACKE_cgetri(LAPACK_ROW_MAJOR, N, &A_inv[0][0], A_inv.cols, &ipiv[0]);
	if (info != 0) throw("LAPACKE_cgetri failed.");

	return A_inv;
}

inline Matrix<complex<double> > Inv(const Matrix<complex<double> > & A_mat)
{
	int N=A_mat.rows;
	if (N!=A_mat.cols) throw("Matrix A_mat must be square.");

	Matrix<complex<double> > A_inv(A_mat);
	valarray<lapack_int> ipiv(N);
	lapack_int info=LAPACKE_zgetrf(LAPACK_ROW_MAJOR, N, N, &A_inv[0][0], A_inv.cols, &ipiv[0]);
	if (info != 0) throw("LAPACKE_zgetrf failed.");
	info=LAPACKE_zgetri(LAPACK_ROW_MAJOR, N, &A_inv[0][0], A_inv.cols, &ipiv[0]);
	if (info != 0) throw("LAPACKE_zgetri failed.");

	return A_inv;
}

inline Matrix<float> SymInv(const Matrix<float> & A_mat)
{
	int N=A_mat.rows;
	if (N!=A_mat.cols) throw("Matrix A_mat must be square.");

	Matrix<float> A_inv(A_mat);
	valarray<lapack_int> ipiv(N);
	lapack_int info=LAPACKE_ssytrf(LAPACK_ROW_MAJOR, 'U', N, &A_inv[0][0], A_inv.cols, &ipiv[0]);
	if (info != 0) throw("LAPACKE_ssytrf failed.");
	info=LAPACKE_ssytri(LAPACK_ROW_MAJOR, 'U', N, &A_inv[0][0], A_inv.cols, &ipiv[0]);
	if (info != 0) throw("LAPACKE_ssytri failed.");

	return A_inv;
}

inline Matrix<double> SymInv(const Matrix<double> & A_mat)
{
	int N=A_mat.rows;
	if (N!=A_mat.cols) throw("Matrix A_mat must be square.");

	Matrix<double> A_inv(A_mat);
	valarray<lapack_int> ipiv(N);
	lapack_int info=LAPACKE_dsytrf(LAPACK_ROW_MAJOR, 'U', N, &A_inv[0][0], A_inv.cols, &ipiv[0]);
	if (info != 0) throw("LAPACKE_dsytrf failed.");
	info=LAPACKE_dsytri(LAPACK_ROW_MAJOR, 'U', N, &A_inv[0][0], A_inv.cols, &ipiv[0]);
	if (info != 0) throw("LAPACKE_dsytri failed.");

	return A_inv;
}

inline Matrix<complex<float> > SymInv(const Matrix<complex<float> > & A_mat)
{
	int N=A_mat.rows;
	if (N!=A_mat.cols) throw("Matrix A_mat must be square.");

	Matrix<complex<float> > A_inv(A_mat);
	valarray<lapack_int> ipiv(N);
	lapack_int info=LAPACKE_csytrf(LAPACK_ROW_MAJOR, 'U', N, &A_inv[0][0], A_inv.cols, &ipiv[0]);
	if (info != 0) throw("LAPACKE_csytrf failed.");
	info=LAPACKE_csytri(LAPACK_ROW_MAJOR, 'U', N, &A_inv[0][0], A_inv.cols, &ipiv[0]);
	if (info != 0) throw("LAPACKE_csytri failed.");

	return A_inv;
}

inline Matrix<complex<double> > SymInv(const Matrix<complex<double> > & A_mat)
{
	int N=A_mat.rows;
	if (N!=A_mat.cols) throw("Matrix A_mat must be square.");

	Matrix<complex<double> > A_inv(A_mat);
	valarray<lapack_int> ipiv(N);
	lapack_int info=LAPACKE_zsytrf(LAPACK_ROW_MAJOR, 'U', N, &A_inv[0][0], A_inv.cols, &ipiv[0]);
	if (info != 0) throw("LAPACKE_zsytrf failed.");
	info=LAPACKE_zsytri(LAPACK_ROW_MAJOR, 'U', N, &A_inv[0][0], A_inv.cols, &ipiv[0]);
	if (info != 0) throw("LAPACKE_zsytri failed.");

	return A_inv;
}

inline Matrix<float> SymPosDefInv(const Matrix<float> & A_mat)
{
	int N=A_mat.rows;
	if (N!=A_mat.cols) throw("Matrix A_mat must be square.");

	Matrix<float> A_inv(A_mat);
	lapack_int info=LAPACKE_spotrf(LAPACK_ROW_MAJOR, 'U', N, &A_inv[0][0], A_mat.cols);
	if (info != 0) throw("LAPACKE_spotrf failed.");
	info=LAPACKE_spotri(LAPACK_ROW_MAJOR, 'U', N, &A_inv[0][0], A_mat.cols);
	if (info != 0) throw("LAPACKE_spotri failed.");

	return A_inv;
}

inline Matrix<double> SymPosDefInv(const Matrix<double> & A_mat)
{
	int N=A_mat.rows;
	if (N!=A_mat.cols) throw("Matrix A_mat must be square.");

	Matrix<double> A_inv(A_mat);
	lapack_int info=LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'U', N, &A_inv[0][0], A_mat.cols);
	if (info != 0) throw("LAPACKE_dpotrf failed.");
	info=LAPACKE_dpotri(LAPACK_ROW_MAJOR, 'U', N, &A_inv[0][0], A_mat.cols);
	if (info != 0) throw("LAPACKE_dpotri failed.");

	return A_inv;
}

inline Matrix<complex<float> > HermPosDefInv(const Matrix<complex<float> > & A_mat)
{
	int N=A_mat.rows;
	if (N!=A_mat.cols) throw("Matrix A_mat must be square.");

	Matrix<complex<float> > A_inv(A_mat);
	lapack_int info=LAPACKE_cpotrf(LAPACK_ROW_MAJOR, 'U', N, &A_inv[0][0], A_mat.cols);
	if (info != 0) throw("LAPACKE_cpotrf failed.");
	info=LAPACKE_cpotri(LAPACK_ROW_MAJOR, 'U', N, &A_inv[0][0], A_mat.cols);
	if (info != 0) throw("LAPACKE_cpotri failed.");

	return A_inv;
}

inline Matrix<complex<double> > HermPosDefInv(const Matrix<complex<double> > & A_mat)
{
	int N=A_mat.rows;
	if (N!=A_mat.cols) throw("Matrix A_mat must be square.");

	Matrix<complex<double> > A_inv(A_mat);
	lapack_int info=LAPACKE_zpotrf(LAPACK_ROW_MAJOR, 'U', N, &A_inv[0][0], A_mat.cols);
	if (info != 0) throw("LAPACKE_zpotrf failed.");
	info=LAPACKE_zpotri(LAPACK_ROW_MAJOR, 'U', N, &A_inv[0][0], A_mat.cols);
	if (info != 0) throw("LAPACKE_zpotri failed.");

	return A_inv;
}

inline void Eig(const Matrix<float> &A_mat, valarray<complex<float> > &EigVal_vec, Matrix<complex<float> > &EigVec_mat)
{
	int N=A_mat.cols;
	if (A_mat.rows!=N) throw("Matrix A_mat must be square.");
	if (EigVal_vec.size()!=N) throw("valarray EigVal_vec must have length same as the number of rows of matrix A_mat.");
	if (EigVec_mat.rows!=N) throw("Matrix EigVec_mat must have the same number of rows as matrix A_mat.");
	if (EigVec_mat.cols!=N) throw("Matrix EigVec_mat must have the same number of columns as matrix A_mat.");
	
	complex<float> i(0,1);
	Matrix<float> A_temp(A_mat), EigVec_r(N,N);
	valarray<float> w_r(N), w_i(N);
	lapack_int info=LAPACKE_sgeev(LAPACK_ROW_MAJOR, 'N', 'V', N, &A_temp[0][0], A_temp.cols, &w_r[0], &w_i[0], &EigVec_r[0][0], EigVec_r.cols, &EigVec_r[0][0], EigVec_r.cols);
	if (info != 0) throw("LAPACKE_sgeev failed.");

	EigVal_vec=w_r+ w_i*i;
	for (int n=0; n<N;n++)
	{
		if (w_i[n]!=0)
		{
			EigVec_mat.col(n)=valarray<float>(EigVec_r.col(n))+i*valarray<float>(EigVec_r.col(n+1));
			EigVec_mat.col(n+1)=valarray<float>(EigVec_r.col(n))-i*valarray<float>(EigVec_r.col(n+1));
			n++;
		}
		else
		{
			EigVec_mat.col(n)=valarray<float>(EigVec_r.col(n))+i*valarray<float>(0.0,N);
		}
	}
}

inline void Eig(const Matrix<double> &A_mat, valarray<complex<double> > &EigVal_vec, Matrix<complex<double> > &EigVec_mat)
{
	int N=A_mat.cols;
	if (A_mat.rows!=N) throw("Matrix A_mat must be square.");
	if (EigVal_vec.size()!=N) throw("valarray EigVal_vec must have length same as the number of rows of matrix A_mat.");
	if (EigVec_mat.rows!=N) throw("Matrix EigVec_mat must have the same number of rows as matrix A_mat.");
	if (EigVec_mat.cols!=N) throw("Matrix EigVec_mat must have the same number of columns as matrix A_mat.");

	complex<double> i(0,1);
	Matrix<double> A_temp(A_mat), EigVec_r(N,N);
	valarray<double> w_r(N), w_i(N);
	lapack_int info=LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', N, &A_temp[0][0], A_temp.cols, &w_r[0], &w_i[0], &EigVec_r[0][0], EigVec_r.cols, &EigVec_r[0][0], EigVec_r.cols);
	if (info != 0) throw("LAPACKE_dgeev failed.");

	EigVal_vec=w_r+ w_i*i;
	for (int n=0; n<N;n++)
	{
		if (w_i[n]!=0)
		{
			EigVec_mat.col(n)=valarray<double>(EigVec_r.col(n))+i*valarray<double>(EigVec_r.col(n+1));
			EigVec_mat.col(n+1)=valarray<double>(EigVec_r.col(n))-i*valarray<double>(EigVec_r.col(n+1));
			n++;
		}
		else
		{
			EigVec_mat.col(n)=valarray<double>(EigVec_r.col(n))+i*valarray<double>(0.0,N);
		}
	}
}

inline void Eig(const Matrix<complex<float> > &A_mat, valarray<complex<float> > &EigVal_vec, Matrix<complex<float> > &EigVec_mat)
{
	int N=A_mat.cols;
	if (A_mat.rows!=N) throw("Matrix A_mat must be square.");
	if (EigVal_vec.size()!=N) throw("valarray EigVal_vec must have length same as the number of rows of matrix A_mat.");
	if (EigVec_mat.rows!=N) throw("Matrix EigVec_mat must have the same number of rows as matrix A_mat.");
	if (EigVec_mat.cols!=N) throw("Matrix EigVec_mat must have the same number of columns as matrix A_mat.");

	complex<float> i(0,1);
	Matrix<complex<float> > A_temp(A_mat);
	EigVal_vec=valarray<complex<float> >(N);
	EigVec_mat=Matrix<complex<float> >(N,N);

	lapack_int info=LAPACKE_cgeev(LAPACK_ROW_MAJOR, 'N', 'V', N, &A_temp[0][0], A_temp.cols, &EigVal_vec[0], &EigVec_mat[0][0], EigVec_mat.cols, &EigVec_mat[0][0], EigVec_mat.cols);
	if (info != 0) throw("LAPACKE_cgeev failed.");
}

inline void Eig(const Matrix<complex<double> > &A_mat, valarray<complex<double> > &EigVal_vec, Matrix<complex<double> > &EigVec_mat)
{
	int N=A_mat.cols;
	if (A_mat.rows!=N) throw("Matrix A_mat must be square.");
	if (EigVal_vec.size()!=N) throw("valarray EigVal_vec must have length same as the number of rows of matrix A_mat.");
	if (EigVec_mat.rows!=N) throw("Matrix EigVec_mat must have the same number of rows as matrix A_mat.");
	if (EigVec_mat.cols!=N) throw("Matrix EigVec_mat must have the same number of columns as matrix A_mat.");

	complex<double> i(0,1);
	Matrix<complex<double> > A_temp(A_mat);
	EigVal_vec=valarray<complex<double> >(N);
	EigVec_mat=Matrix<complex<double> >(N,N);

	lapack_int info=LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'N', 'V', N, &A_temp[0][0], A_temp.cols, &EigVal_vec[0], &EigVec_mat[0][0], EigVec_mat.cols, &EigVec_mat[0][0], EigVec_mat.cols);
	if (info != 0) throw("LAPACKE_zgeev failed.");
}

inline void SymEig(const Matrix<float> &A_mat, valarray<float> &EigVal_vec, Matrix<float> &EigVec_mat)
{
	int N=A_mat.cols;
	if (A_mat.rows!=N) throw("Matrix A_mat must be square.");
	if (EigVal_vec.size()!=N) throw("valarray EigVal_vec must have length same as the number of rows of matrix A_mat.");
	if (EigVec_mat.rows!=N) throw("Matrix EigVec_mat must have the same number of rows as matrix A_mat.");
	if (EigVec_mat.cols!=N) throw("Matrix EigVec_mat must have the same number of columns as matrix A_mat.");

	EigVec_mat=A_mat;
	EigVal_vec=valarray<float>(N);

	lapack_int info=LAPACKE_ssyevd(LAPACK_ROW_MAJOR, 'V', 'U', N, &EigVec_mat[0][0], EigVec_mat.cols, &EigVal_vec[0]);
	if (info != 0) throw("LAPACKE_ssyevd failed.");
}

inline void SymEig(const Matrix<double> &A_mat, valarray<double> &EigVal_vec, Matrix<double> &EigVec_mat)
{
	int N=A_mat.cols;
	if (A_mat.rows!=N) throw("Matrix A_mat must be square.");
	if (EigVal_vec.size()!=N) throw("valarray EigVal_vec must have length same as the number of rows of matrix A_mat.");
	if (EigVec_mat.rows!=N) throw("Matrix EigVec_mat must have the same number of rows as matrix A_mat.");
	if (EigVec_mat.cols!=N) throw("Matrix EigVec_mat must have the same number of columns as matrix A_mat.");

	EigVec_mat=A_mat;
	EigVal_vec=valarray<double>(N);

	lapack_int info=LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'V', 'U', N, &EigVec_mat[0][0], EigVec_mat.cols, &EigVal_vec[0]);
	if (info != 0) throw("LAPACKE_dsyevd failed.");
}

inline void HermEig(const Matrix<complex<float> > &A_mat, valarray<float> &EigVal_vec, Matrix<complex<float> > &EigVec_mat)
{
	int N=A_mat.cols;
	if (A_mat.rows!=N) throw("Matrix A_mat must be square.");
	if (EigVal_vec.size()!=N) throw("valarray EigVal_vec must have length same as the number of rows of matrix A_mat.");
	if (EigVec_mat.rows!=N) throw("Matrix EigVec_mat must have the same number of rows as matrix A_mat.");
	if (EigVec_mat.cols!=N) throw("Matrix EigVec_mat must have the same number of columns as matrix A_mat.");

	EigVec_mat=A_mat;
	EigVal_vec=valarray<float>(N);

	lapack_int info=LAPACKE_cheevd(LAPACK_ROW_MAJOR, 'V', 'U', N, &EigVec_mat[0][0], EigVec_mat.cols, &EigVal_vec[0]);
	if (info != 0) throw("LAPACKE_cheevd failed.");
}

inline void HermEig(const Matrix<complex<double> > &A_mat, valarray<double> &EigVal_vec, Matrix<complex<double> > &EigVec_mat)
{
	int N=A_mat.cols;
	if (A_mat.rows!=N) throw("Matrix A_mat must be square.");
	if (EigVal_vec.size()!=N) throw("valarray EigVal_vec must have length same as the number of rows of matrix A_mat.");
	if (EigVec_mat.rows!=N) throw("Matrix EigVec_mat must have the same number of rows as matrix A_mat.");
	if (EigVec_mat.cols!=N) throw("Matrix EigVec_mat must have the same number of columns as matrix A_mat.");

	EigVec_mat=A_mat;
	EigVal_vec=valarray<double>(N);

	lapack_int info=LAPACKE_zheevd(LAPACK_ROW_MAJOR, 'V', 'U', N, &EigVec_mat[0][0], EigVec_mat.cols, &EigVal_vec[0]);
	if (info != 0) throw("LAPACKE_zheevd failed.");
}

inline void GenEig(const Matrix<float> &A_mat, const Matrix<float> &B_mat, valarray<complex<float> > &EigVal_vec, Matrix<complex<float> > &EigVec_mat)
{
	int N=A_mat.cols;
	if (A_mat.rows!=N) throw("Matrix A_mat must be square.");
	if (B_mat.cols!=B_mat.rows) throw("Matrix B_mat must be square.");
	if (N!=B_mat.cols) throw("Matrix A_mat and B_mat must have the same size.");
	if (EigVal_vec.size()!=N) throw("valarray EigVal_vec must have length same as the number of rows of matrix A_mat.");
	if (EigVec_mat.rows!=N) throw("Matrix EigVec_mat must have the same number of rows as matrix A_mat.");
	if (EigVec_mat.cols!=N) throw("Matrix EigVec_mat must have the same number of columns as matrix A_mat.");
	
	complex<float> i(0,1);
	Matrix<float> A_temp(A_mat), B_temp(B_mat), EigVec_r(N,N);
	valarray<float> alphar(N), alphai(N), beta(N);

	lapack_int info=LAPACKE_sggev(LAPACK_ROW_MAJOR, 'N', 'V', N, &A_temp[0][0], A_temp.cols, &B_temp[0][0], B_temp.cols, &alphar[0], &alphai[0], &beta[0], &EigVec_r[0][0], EigVec_r.cols, &EigVec_r[0][0], EigVec_r.cols);
	if (info != 0) throw("LAPACKE_sggev failed.");

	EigVal_vec=(alphar+ alphai*i)/beta;
	for (int n=0; n<N;n++)
	{
		if (alphai[n]!=0)
		{
			EigVec_mat.col(n)=valarray<float>(EigVec_r.col(n))+i*valarray<float>(EigVec_r.col(n+1));
			EigVec_mat.col(n+1)=valarray<float>(EigVec_r.col(n))-i*valarray<float>(EigVec_r.col(n+1));
			n++;
		}
		else
		{
			EigVec_mat.col(n)=valarray<float>(EigVec_r.col(n))+i*valarray<float>(0.0,N);
		}
	}
}

inline void GenEig(const Matrix<double> &A_mat, const Matrix<double> &B_mat, valarray<complex<double> > &EigVal_vec, Matrix<complex<double> > &EigVec_mat)
{
	int N=A_mat.cols;
	if (A_mat.rows!=N) throw("Matrix A_mat must be square.");
	if (B_mat.cols!=B_mat.rows) throw("Matrix B_mat must be square.");
	if (N!=B_mat.cols) throw("Matrix A_mat and B_mat must have the same size.");
	if (EigVal_vec.size()!=N) throw("valarray EigVal_vec must have length same as the number of rows of matrix A_mat.");
	if (EigVec_mat.rows!=N) throw("Matrix EigVec_mat must have the same number of rows as matrix A_mat.");
	if (EigVec_mat.cols!=N) throw("Matrix EigVec_mat must have the same number of columns as matrix A_mat.");

	complex<double> i(0,1);
	Matrix<double> A_temp(A_mat), B_temp(B_mat), EigVec_r(N,N);
	valarray<double> alphar(N), alphai(N), beta(N);

	lapack_int info=LAPACKE_dggev(LAPACK_ROW_MAJOR, 'N', 'V', N, &A_temp[0][0], A_temp.cols, &B_temp[0][0], B_temp.cols, &alphar[0], &alphai[0], &beta[0], &EigVec_r[0][0], EigVec_r.cols, &EigVec_r[0][0], EigVec_r.cols);
	if (info != 0) throw("LAPACKE_dggev failed.");

	EigVal_vec=(alphar+alphai*i)/beta;
	for (int n=0; n<N;n++)
	{
		if (alphai[n]!=0)
		{
			EigVec_mat.col(n)=valarray<double>(EigVec_r.col(n))+i*valarray<double>(EigVec_r.col(n+1));
			EigVec_mat.col(n+1)=valarray<double>(EigVec_r.col(n))-i*valarray<double>(EigVec_r.col(n+1));
			n++;
		}
		else
		{
			EigVec_mat.col(n)=valarray<double>(EigVec_r.col(n))+i*valarray<double>(0.0,N);
		}
	}
}

inline void GenEig(const Matrix<complex<float> > &A_mat, const Matrix<complex<float> > &B_mat, valarray<complex<float> > &EigVal_vec, Matrix<complex<float> > &EigVec_mat)
{
	int N=A_mat.cols;
	if (A_mat.rows!=N) throw("Matrix A_mat must be square.");
	if (B_mat.cols!=B_mat.rows) throw("Matrix B_mat must be square.");
	if (N!=B_mat.cols) throw("Matrix A_mat and B_mat must have the same size.");
	if (EigVal_vec.size()!=N) throw("valarray EigVal_vec must have length same as the number of rows of matrix A_mat.");
	if (EigVec_mat.rows!=N) throw("Matrix EigVec_mat must have the same number of rows as matrix A_mat.");
	if (EigVec_mat.cols!=N) throw("Matrix EigVec_mat must have the same number of columns as matrix A_mat.");

	complex<float> i(0,1);
	Matrix<complex<float> > A_temp(A_mat), B_temp(B_mat);
	valarray<complex<float> > alpha(N), beta(N);

	lapack_int info=LAPACKE_cggev(LAPACK_ROW_MAJOR, 'N', 'V', N, &A_temp[0][0], A_temp.cols, &B_temp[0][0], B_temp.cols, &alpha[0], &beta[0], &EigVec_mat[0][0], EigVec_mat.cols, &EigVec_mat[0][0], EigVec_mat.cols);
	if (info != 0) throw("LAPACKE_cggev failed.");

	EigVal_vec=alpha/beta;
}

inline void GenEig(const Matrix<complex<double> > &A_mat, const Matrix<complex<double> > &B_mat, valarray<complex<double> > &EigVal_vec, Matrix<complex<double> > &EigVec_mat)
{
	int N=A_mat.cols;
	if (A_mat.rows!=N) throw("Matrix A_mat must be square.");
	if (B_mat.cols!=B_mat.rows) throw("Matrix B_mat must be square.");
	if (N!=B_mat.cols) throw("Matrix A_mat and B_mat must have the same size.");
	if (EigVal_vec.size()!=N) throw("valarray EigVal_vec must have length same as the number of rows of matrix A_mat.");
	if (EigVec_mat.rows!=N) throw("Matrix EigVec_mat must have the same number of rows as matrix A_mat.");
	if (EigVec_mat.cols!=N) throw("Matrix EigVec_mat must have the same number of columns as matrix A_mat.");

	complex<double> i(0,1);
	Matrix<complex<double> > A_temp(A_mat), B_temp(B_mat);
	valarray<complex<double> > alpha(N), beta(N);

	lapack_int info=LAPACKE_zggev(LAPACK_ROW_MAJOR, 'N', 'V', N, &A_temp[0][0], A_temp.cols, &B_temp[0][0], B_temp.cols, &alpha[0], &beta[0], &EigVec_mat[0][0], EigVec_mat.cols, &EigVec_mat[0][0], EigVec_mat.cols);
	if (info != 0) throw("LAPACKE_zggev failed.");

	EigVal_vec=alpha/beta;
}

inline void GenSymDefEig(const Matrix<float> &A_mat, const Matrix<float> &B_mat, valarray<float> &EigVal_vec, Matrix<float> &EigVec_mat)
{
	int N=A_mat.cols;
	if (A_mat.rows!=N) throw("Matrix A_mat must be square.");
	if (B_mat.cols!=B_mat.rows) throw("Matrix B_mat must be square.");
	if (N!=B_mat.cols) throw("Matrix A_mat and B_mat must have the same size.");
	if (EigVal_vec.size()!=N) throw("valarray EigVal_vec must have length same as the number of rows of matrix A_mat.");
	if (EigVec_mat.rows!=N) throw("Matrix EigVec_mat must have the same number of rows as matrix A_mat.");
	if (EigVec_mat.cols!=N) throw("Matrix EigVec_mat must have the same number of columns as matrix A_mat.");

	Matrix<float> B_temp(B_mat);

	EigVal_vec=valarray<float>(N);
	EigVec_mat=A_mat;
	lapack_int info=LAPACKE_ssygvd(LAPACK_ROW_MAJOR, 1, 'V', 'U', N, &EigVec_mat[0][0], EigVec_mat.cols, &B_temp[0][0], B_temp.cols, &EigVal_vec[0]);
	if (info > N)
		throw("Error using LAPACKE_ssygvd, matrix B_mat is not positive definite.");
	else if (info != 0)
		throw("LAPACKE_ssygvd failed.");
}

inline void GenSymDefEig(const Matrix<double> &A_mat, const Matrix<double> &B_mat, valarray<double> &EigVal_vec, Matrix<double> &EigVec_mat)
{
	int N=A_mat.cols;
	if (A_mat.rows!=N) throw("Matrix A_mat must be square.");
	if (B_mat.cols!=B_mat.rows) throw("Matrix B_mat must be square.");
	if (N!=B_mat.cols) throw("Matrix A_mat and B_mat must have the same size.");
	if (EigVal_vec.size()!=N) throw("valarray EigVal_vec must have length same as the number of rows of matrix A_mat.");
	if (EigVec_mat.rows!=N) throw("Matrix EigVec_mat must have the same number of rows as matrix A_mat.");
	if (EigVec_mat.cols!=N) throw("Matrix EigVec_mat must have the same number of columns as matrix A_mat.");

	Matrix<double> B_temp(B_mat);

	EigVal_vec=valarray<double>(N);
	EigVec_mat=A_mat;
	lapack_int info=LAPACKE_dsygvd(LAPACK_ROW_MAJOR, 1, 'V', 'U', N, &EigVec_mat[0][0], EigVec_mat.cols, &B_temp[0][0], B_temp.cols, &EigVal_vec[0]);
	if (info > N)
		throw("Error using LAPACKE_dsygvd, matrix B_mat is not positive definite.");
	else if (info != 0)
		throw("LAPACKE_dsygvd failed.");
}

inline void GenHermDefEig(const Matrix<complex<float> > &A_mat, const Matrix<complex<float> > &B_mat, valarray<float> &EigVal_vec, Matrix<complex<float> > &EigVec_mat)
{
	int N=A_mat.cols;
	if (A_mat.rows!=N) throw("Matrix A_mat must be square.");
	if (B_mat.cols!=B_mat.rows) throw("Matrix B_mat must be square.");
	if (N!=B_mat.cols) throw("Matrix A_mat and B_mat must have the same size.");
	if (EigVal_vec.size()!=N) throw("valarray EigVal_vec must have length same as the number of rows of matrix A_mat.");
	if (EigVec_mat.rows!=N) throw("Matrix EigVec_mat must have the same number of rows as matrix A_mat.");
	if (EigVec_mat.cols!=N) throw("Matrix EigVec_mat must have the same number of columns as matrix A_mat.");

	Matrix<complex<float> > B_temp(B_mat);

	EigVal_vec=valarray<float>(N);
	EigVec_mat=A_mat;
	lapack_int info=LAPACKE_chegvd(LAPACK_ROW_MAJOR, 1, 'V', 'U', N, &EigVec_mat[0][0], EigVec_mat.cols, &B_temp[0][0], B_temp.cols, &EigVal_vec[0]);
	if (info > N)
		throw("Error using LAPACKE_chegvd, matrix B_mat is not positive definite.");
	else if (info != 0)
		throw("LAPACKE_chegvd failed.");
}

inline void GenHermDefEig(const Matrix<complex<double> > &A_mat, const Matrix<complex<double> > &B_mat, valarray<double> &EigVal_vec, Matrix<complex<double> > &EigVec_mat)
{
	int N=A_mat.cols;
	if (A_mat.rows!=N) throw("Matrix A_mat must be square.");
	if (B_mat.cols!=B_mat.rows) throw("Matrix B_mat must be square.");
	if (N!=B_mat.cols) throw("Matrix A_mat and B_mat must have the same size.");
	if (EigVal_vec.size()!=N) throw("valarray EigVal_vec must have length same as the number of rows of matrix A_mat.");
	if (EigVec_mat.rows!=N) throw("Matrix EigVec_mat must have the same number of rows as matrix A_mat.");
	if (EigVec_mat.cols!=N) throw("Matrix EigVec_mat must have the same number of columns as matrix A_mat.");

	Matrix<complex<double> > B_temp(B_mat);

	EigVal_vec=valarray<double>(N);
	EigVec_mat=A_mat;
	lapack_int info=LAPACKE_zhegvd(LAPACK_ROW_MAJOR, 1, 'V', 'U', N, &EigVec_mat[0][0], EigVec_mat.cols, &B_temp[0][0], B_temp.cols, &EigVal_vec[0]);
	if (info > N)
		throw("Error using LAPACKE_zhegvd, matrix B_mat is not positive definite.");
	else if (info != 0)
		throw("LAPACKE_zhegvd failed.");
}

#endif