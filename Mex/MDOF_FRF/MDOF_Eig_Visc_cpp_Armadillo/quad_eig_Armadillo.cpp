#include <complex>
#include <math.h>

#include "ArmadilloMatrixOperations.hpp"

using namespace std;
using namespace arma;

void quad_eig(const Mat<double>& K_mat,const Mat<double>& C_mat,const Mat<double>& M_mat, Mat<complex<double> >& Epsi_normalized_mat, Col<complex<double> >& EigValues_col)
{
	size_t N=M_mat.n_rows;
	Mat<double> N_mat=M_mat;
	Mat<double> A_mat=join_vert(join_horiz(-K_mat,zeros<Mat<double > >(N,N)),join_horiz(zeros<Mat<double > >(N,N),N_mat));
	Mat<double> B_mat=join_vert(join_horiz(C_mat,M_mat),join_horiz(N_mat, zeros<Mat<double > >(N,N)));

	Mat<complex<double> > Phi_mat;
	eig_pair(EigValues_col, Phi_mat, A_mat, B_mat);

	Col<complex<double> > B_r_vec=diagvec((strans(Phi_mat)*B_mat*Phi_mat));
	Mat<complex<double> > Phi_normalized_mat(2*N,2*N);
	for (size_t n=0;n<2*N;n++)
		Phi_normalized_mat.col(n)=Phi_mat.col(n)/sqrt(B_r_vec(n));

	Epsi_normalized_mat=Phi_normalized_mat.rows(0,N-1);
}