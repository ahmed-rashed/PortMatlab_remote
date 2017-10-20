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

void quad_eig(const Matrix<double>& K_mat,const Matrix<double>& C_mat,const Matrix<double>& M_mat, Matrix<complex<double> >& Epsi_normalized_mat, valarray<complex<double> >& EigVal_vec)
{
	const complex<double> i(0,1);
	size_t N=M_mat.cols;
	if (EigVal_vec.size()!=2*N) throw("valarray EigVal_vec must have length same as the number of rows of matrix M_mat.");
	if (Epsi_normalized_mat.rows!=N) throw("Matrix Epsi_normalized_mat must have twice the number of rows as matrix M_mat.");
	if (Epsi_normalized_mat.cols!=2*N) throw("Matrix Epsi_normalized_mat must have twice the number of columns as matrix M_mat.");

	Matrix<double> N_mat(M_mat);
	
	Matrix<double> A_mat(2*N,2*N,0.0);
	A_mat.SubMatrix(0,N-1,0,N-1)=valarray<double>(-K_mat.SubMatrix(0,N-1,0,N-1));
	A_mat.SubMatrix(N,2*N-1,N,2*N-1)=valarray<double>(N_mat.SubMatrix(0,N-1,0,N-1));

	Matrix<double> B_mat(2*N,2*N,0.0);
	B_mat.SubMatrix(0,N-1,0,N-1)=valarray<double>(C_mat.SubMatrix(0,N-1,0,N-1));
	B_mat.SubMatrix(0,N-1,N,2*N-1)=valarray<double>(M_mat.SubMatrix(0,N-1,0,N-1));
	B_mat.SubMatrix(N,2*N-1,0,N-1)=valarray<double>(N_mat.SubMatrix(0,N-1,0,N-1));

	Matrix<complex<double> > Phi_mat(2*N,2*N), Phi_normalized_mat(2*N,2*N);
	GenEig(A_mat, B_mat, EigVal_vec, Phi_mat);
	valarray<complex<double> > B_r_vec(MatDiag(MatMul(MatMul(transpose(Phi_mat),B_mat),Phi_mat)));

	Phi_normalized_mat=scale_cols(Phi_mat,1.0/sqrt(B_r_vec));

	Epsi_normalized_mat.SubMatrix(0,N-1,0,2*N-1)=valarray<complex<double> >(Phi_normalized_mat.SubMatrix(0,N-1,0,2*N-1));
}