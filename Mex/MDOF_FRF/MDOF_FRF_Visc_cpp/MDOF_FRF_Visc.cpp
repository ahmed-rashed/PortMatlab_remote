#include <complex>
#include <valarray>
#include <math.h>
#include "valarray_ComplexRealOperators.hpp"
#include "Matrix_BLAS.hpp"
#include "MatrixOperations.hpp"
#include "a_matrices.hpp"

using namespace std;
using namespace a_matrices;

Matrix<complex<double> > MDOF_FRF_Visc(const valarray<complex<double> >& EigValues_vec, const Matrix<complex<double> >& EigVectors_Normalized_mat, const Matrix<double>& w_column, const valarray<size_t>& m_vec, const valarray<size_t>& n_vec)
{
	size_t N=EigVectors_Normalized_mat.rows;
	size_t N_w=w_column.size();
	size_t N_cols=m_vec.size();
	const complex<double> i(0,1);

	Matrix<complex<double> > H_w_mat_SDOF(N_w,N_cols,complex<double>(0,0));
	Matrix<complex<double> > H_w_n_m_cols(N_w,N_cols,complex<double>(0,0));
	valarray<size_t> A_ind=RowMajorSub2Ind(N, m_vec, n_vec);

	for (int r=0;r<2*N;r++)
	{
		Matrix<complex<double> > EigVector_r(N,1,valarray<complex<double> >(EigVectors_Normalized_mat.col(r)));
		Matrix<complex<double> > A_r_mat=MatMul(EigVector_r,transpose(EigVector_r));
		Matrix<complex<double> > A_r_temp_row(1,N_cols,valarray<complex<double> >(A_r_mat[A_ind]));
		H_w_mat_SDOF=H_w_mat_SDOF+MatMul(1.0/(i*w_column-EigValues_vec[r]),A_r_temp_row);

		if ((imag(EigValues_vec[r])!=0) && ((r+1)%2!=0))
			continue;

		H_w_n_m_cols=H_w_n_m_cols+H_w_mat_SDOF;
		H_w_mat_SDOF= complex<double>(0, 0);
	}

	return H_w_n_m_cols;
}

