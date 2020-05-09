#include <complex>
#include <math.h>

#include "ArmadilloMatrixOperations.hpp"

using namespace std;
using namespace arma;

void quad_eig(const Mat<double>& K_mat,const Mat<double>& C_mat,const Mat<double>& M_mat, Mat<complex<double> >& Epsi_normalized_mat, Col<complex<double> >& EigValues_col);

void MDOF_Eig_Visc(const Mat<double>& M_mat,const Mat<double>& C_mat,const Mat<double>& K_mat, const bool& isPropotional, Mat<complex<double> >& EigVectors_Normalized_mat, Col<complex<double> >& EigValues_col)
{
	const complex<double> i(0,1);
	size_t	N=M_mat.n_rows;
	uvec nn=conv_to<uvec>::from(stepvec(1,2,2*N-1)-1);
	if(isPropotional || all(all(C_mat==0)))    //Undamped or proportional
	{
		Col<complex<double> > EigValues_U_vec;
		Mat<complex<double> > EigVectors_U_mat;
		eig_pair(EigValues_U_vec, EigVectors_U_mat, -K_mat, M_mat);

		Col<double> M_r_col=diagvec(strans(real(EigVectors_U_mat))*M_mat*real(EigVectors_U_mat));
		Col<double> w_U_r_col=sqrt(-real(EigValues_U_vec));

		Col<double> C_r_col=diagvec(strans(real(EigVectors_U_mat))*C_mat*real(EigVectors_U_mat));
		Col<double> zeta_col=C_r_col/(2*M_r_col%w_U_r_col);
		Col<double> w_d_r_col=w_U_r_col%sqrt(1-square(zeta_col));
		Col<complex<double> > EigValues_col_temp=-zeta_col%w_U_r_col+i*w_d_r_col;

		EigValues_col=Col<complex<double> >(2*N);
		EigValues_col(nn)=EigValues_col_temp;
		EigValues_col(nn+1)=conj(EigValues_col_temp);

		EigVectors_Normalized_mat=Mat<complex<double> >(N,2*N);
		EigVectors_Normalized_mat.cols(nn)=real(EigVectors_U_mat)*diagmat(1/(sqrt(i*2.0*w_d_r_col%M_r_col)));
		EigVectors_Normalized_mat.cols(nn+1)=conj(EigVectors_Normalized_mat.cols(nn));
	}
	else 
	{
		quad_eig(K_mat, C_mat, M_mat, EigVectors_Normalized_mat, EigValues_col);
	}

	uvec Index=sort_index(abs(imag(EigValues_col(nn))));
	uvec IIndex(2*N);
	IIndex(nn)=2*Index;
	IIndex(nn+1)=2*Index+1;
	EigValues_col=EigValues_col(IIndex);
	EigVectors_Normalized_mat=EigVectors_Normalized_mat.cols(IIndex);

	return;
}
