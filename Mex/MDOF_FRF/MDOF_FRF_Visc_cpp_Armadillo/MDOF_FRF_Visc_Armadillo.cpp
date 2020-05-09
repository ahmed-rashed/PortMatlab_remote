#include <complex>

#include <math.h>

#include "ArmadilloMatrixOperations.hpp"
//#include "armaMex.hpp"

using namespace std;
using namespace arma;

Mat<complex<double> > MDOF_FRF_Visc(const Col<complex<double> >& EigValues_col, const Mat<complex<double> >& EigVectors_Normalized_mat, const Mat<double>& w_column, const Row<double>& n_C, const Row<double>& m_C)
{
	const size_t N=EigVectors_Normalized_mat.n_rows;
	const size_t n_w=w_column.size();

	//The following check was already performed in the mexFunction. Hence it is unnecessary here
	//if (n_C.size() != m_C.size()) throw("Dimensions of n_C and m_C must be identical");

	const size_t n_col=n_C.size();
	const complex<double> i(0,1);

	Mat<complex<double> > H_w_mat_SDOF(n_w,n_col);
	H_w_mat_SDOF.zeros();
	Mat<complex<double> > H_w_n_m_cols(n_w,n_col);
	H_w_n_m_cols.zeros();
	Row<uword> A_ind_row=ColMajorSub2Ind(N, n_C, m_C);	//Armadillo is column major
	for (int r=0;r<2*N;r++)
	{
		Mat<complex<double> > A_r_mat=EigVectors_Normalized_mat.col(r)*strans(EigVectors_Normalized_mat.col(r));
		Row<complex<double> > A_r_temp_row=strans(A_r_mat(A_ind_row));
		H_w_mat_SDOF=H_w_mat_SDOF+(1/(i*w_column-EigValues_col(r)))*A_r_temp_row;

		if ((imag(EigValues_col(r))!=0) && ((r+1)%2!=0))
			continue;

		H_w_n_m_cols=H_w_n_m_cols+H_w_mat_SDOF;
		H_w_mat_SDOF.zeros();
	}

	return H_w_n_m_cols;
}

