#include <complex>
#include <valarray>
#include <math.h>
#include <algorithm>
#include <numeric>
#include "valarray_ComplexRealOperators.h"
#include "Matrix_BLAS.h"
#include "Matrix_LAPACK.h"
#include "a_matrices.h"

using namespace std;
using namespace a_matrices;

void quad_eig(const Matrix<double>& K_mat,const Matrix<double>& C_mat,const Matrix<double>& M_mat, Matrix<complex<double> >& Epsi_normalized_mat, valarray<complex<double> >& EigVal_vec);

void MDOF_Eig_Visc(const Matrix<double>& M_mat,const Matrix<double>& C_mat,const Matrix<double>& K_mat, const bool& isPropotional, Matrix<complex<double> >& EigVectors_Normalized, valarray<complex<double> >& EigValues_vec)
{
	size_t N=M_mat.cols;
	if (EigValues_vec.size()!=2*N) throw("valarray EigValues_vec must have length same as the number of rows of matrix M_mat.");
	if (EigVectors_Normalized.rows!=N) throw("Matrix Epsi_normalized_mat must have the same number of rows as matrix M_mat.");
	if (EigVectors_Normalized.cols!=2*N) throw("Matrix Epsi_normalized_mat must have twice the number of columns as matrix M_mat.");

	const complex<double> i(0,1);
	if (isPropotional || all(C_mat==0.0))    //Undamped or proportional
	{
		Matrix<double> EigVectors_U_mat(N,N);
		valarray<double> EigValues_U_vec(N);
		GenSymDefEig(-K_mat,M_mat, EigValues_U_vec, EigVectors_U_mat);
    
		valarray<double> M_r_vec(MatDiag(MatMul(MatMul(transpose(EigVectors_U_mat),M_mat),EigVectors_U_mat)));
		valarray<double> w_U_r_vec(sqrt(-EigValues_U_vec));

		valarray<double> C_r_vec=MatDiag(MatMul(MatMul(transpose(EigVectors_U_mat),C_mat),EigVectors_U_mat));
		valarray<double> zeta_r_vec(C_r_vec/(2.0*M_r_vec*w_U_r_vec));
		valarray<double> w_d_r_vec(w_U_r_vec*sqrt(1.0-zeta_r_vec*zeta_r_vec));
		valarray<complex<double> > EigValues_vec_temp(-zeta_r_vec*w_U_r_vec+i*w_d_r_vec);

		EigValues_vec[slice(0,N,2)]=EigValues_vec_temp;
		EigValues_vec[slice(1,N,2)]=conj(EigValues_vec_temp);
    
		EigVectors_Normalized.SubMatrix(0, N-1, 1, 0, 2*N-2, 2)=valarray<complex<double> >(scale_cols(EigVectors_U_mat+i*Matrix<double>(N,N,0.0),1.0/sqrt(i*2.0*w_d_r_vec*M_r_vec)).SubMatrix(0, N-1, 0, N-1));
		EigVectors_Normalized.SubMatrix(0, N-1, 1, 1, 2*N-1, 2)=conj(valarray<complex<double> >(EigVectors_Normalized.SubMatrix(0, N-1, 1, 0, 2*N-2, 2)));
	}
	else    //Non-proportional
	{
		quad_eig(K_mat,C_mat,M_mat, EigVectors_Normalized, EigValues_vec);
	}

	valarray<size_t> Index(N);
	iota(begin(Index), end(Index), 0);
	Index=size_t(2)*Index;
	sort(begin(Index), end(Index), [&](size_t a, size_t b) { return abs(imag(EigValues_vec[a])) < abs(imag(EigValues_vec[b])); });
	valarray<size_t> IIndex(2*N);
	IIndex[slice(0,N,2)]=Index;
	IIndex[slice(1,N,2)]=Index+size_t(1);
	valarray<complex<double> > EigValuesTemp_vec(EigValues_vec);
	EigValues_vec=EigValuesTemp_vec[IIndex];
	Matrix<complex<double> > EigVectors_NormalizedTemp(EigVectors_Normalized);
	for (size_t n=0;n<2*N;n++)
		EigVectors_Normalized.col(n)=valarray<complex<double> >(EigVectors_NormalizedTemp.col(IIndex[n]));

	return;
}