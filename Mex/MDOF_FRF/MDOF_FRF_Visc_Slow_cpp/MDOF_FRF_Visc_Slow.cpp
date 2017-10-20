#include <complex>
#include <valarray>
#include <math.h>
#include "valarray_ComplexRealOperators.h"
#include "Matrix_LAPACK.h"
#include "a_matrices.h"

using namespace std;
using namespace a_matrices;

Matrix<complex<double> > MDOF_FRF_Visc_Slow(const Matrix<double>& M_mat,const Matrix<double>& C_mat,const Matrix<double>& K_mat,const Matrix<double>& w_col, const valarray<size_t>& index_vec)
{
	const double pi=atan(1.0)*4;
	const complex<double> i(0,1);
	const size_t n_w_points=w_col.size();
	const size_t n_index_vec=index_vec.size();

	Matrix<complex<double> > H_cols(n_w_points,n_index_vec);

	for (unsigned int n=0; n<n_w_points; n++)
	{
		Matrix<complex<double> > H_w_mat=SymInv(-w_col[n][0]*w_col[n][0]*M_mat+i*w_col[n][0]*C_mat+K_mat);

		H_cols.row(n)=H_w_mat[index_vec];
	}

	return H_cols;
}
