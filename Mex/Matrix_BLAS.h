#ifndef Matrix_BLAS_h
#define Matrix_BLAS_h

#include <valarray>
#include <complex>
#define MKL_Complex8 std::complex<float>
#define MKL_Complex16 std::complex<double>
#include "a_matrices.h"
#include "mkl_cblas.h"

using namespace std;
using namespace a_matrices;

inline Matrix<float> MatMul(const Matrix<float> &lhs, const Matrix<float> &rhs)
{
	if (lhs.cols!=rhs.rows) throw("Inner matrices dimensions must match.");
	Matrix<float> Res(lhs.rows,rhs.cols,0.0);
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, lhs.rows, rhs.cols, lhs.cols, 1.0, &lhs[0][0], lhs.cols, &rhs[0][0], rhs.cols, 0.0, &Res[0][0], Res.cols);
	return Res;
}

inline Matrix<double> MatMul(const Matrix<double> &lhs, const Matrix<double> &rhs)
{
	if (lhs.cols!=rhs.rows) throw("Inner matrices dimensions must match.");
	Matrix<double> Res(lhs.rows,rhs.cols,0.0);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, lhs.rows, rhs.cols, lhs.cols, 1.0, &lhs[0][0], lhs.cols, &rhs[0][0], rhs.cols, 0.0, &Res[0][0], Res.cols);
	return Res;
}

inline Matrix<complex<float> > MatMul(const Matrix<complex<float> > &lhs, const Matrix<complex<float> > &rhs)
{
	if (lhs.cols!=rhs.rows) throw("Inner matrices dimensions must match.");
	const complex<float> zero(0.0,0.0);
	const complex<float> one(1.0,0.0);
	Matrix<complex<float> > Res(lhs.rows,rhs.cols,zero);
	cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, lhs.rows, rhs.cols, lhs.cols, &one, &lhs[0][0], lhs.cols, &rhs[0][0], rhs.cols, &zero, &Res[0][0], Res.cols);
	return Res;
}

inline Matrix<complex<double> > MatMul(const Matrix<complex<double> > &lhs, const Matrix<complex<double> > &rhs)
{
	if (lhs.cols!=rhs.rows) throw("Inner matrices dimensions must match.");
	const complex<double> zero(0.0,0.0);
	const complex<double> one(1.0,0.0);
	Matrix<complex<double> > Res(lhs.rows,rhs.cols,zero);
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, lhs.rows, rhs.cols, lhs.cols, &one, &lhs[0][0], lhs.cols, &rhs[0][0], rhs.cols, &zero, &Res[0][0], Res.cols);
	return Res;
}

inline Matrix<complex<float> > MatMul(const Matrix<float> &lhs, const Matrix<complex<float> > &rhs)
{
	const complex<float> i(0,1);
	Matrix<complex<float> > llhs(lhs+i*Matrix<float>(lhs.rows,lhs.cols,0.0));
	return MatMul(llhs, rhs);
}

inline Matrix<complex<float> > MatMul(const Matrix<complex<float> > &lhs, const Matrix<float> &rhs)
{
	const complex<float> i(0,1);
	Matrix<complex<float> > rrhs(rhs+i*Matrix<float>(rhs.rows,rhs.cols,0.0));
	return MatMul(lhs, rrhs);
}

inline Matrix<complex<double> > MatMul(const Matrix<double> &lhs, const Matrix<complex<double> > &rhs)
{
	const complex<double> i(0,1);
	Matrix<complex<double> > llhs(lhs+i*Matrix<double>(lhs.rows,lhs.cols,0.0));
	return MatMul(llhs, rhs);
}

inline Matrix<complex<double> > MatMul(const Matrix<complex<double> > &lhs, const Matrix<double> &rhs)
{
	const complex<double> i(0,1);
	Matrix<complex<double> > rrhs(rhs+i*Matrix<double>(rhs.rows,rhs.cols,0.0));
	return MatMul(lhs, rrhs);
}

#endif