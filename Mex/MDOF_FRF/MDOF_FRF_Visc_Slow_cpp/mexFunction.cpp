#include <complex>
#include <valarray>
#include <math.h>
#include "mex.h"
//#include "ComplexRealOperators.h"
#include "MatrixOperations.h"
#include "a_matrices.h"

using namespace std;
using namespace a_matrices;

Matrix<complex<double> > MDOF_FRF_Visc_Slow(const Matrix<double>& M_mat,const Matrix<double>& C_mat,const Matrix<double>& K_mat,const Matrix<double>& w_col, const valarray<size_t>& index_vec);

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 
    // Check for proper number of arguments
    if (nrhs != 6)
		mexErrMsgTxt("Dear student, This function needs 6 inputs."); 

	if (nlhs > 1)
		mexErrMsgTxt("Dear student, this function returns only one output."); 
    
	size_t N=mxGetM(prhs[0]);
	if (N != mxGetN(prhs[0]))
		mexErrMsgTxt("Dear student, The M_mat matrix must be square.");

	for (unsigned int n=1; n<3; n++)
	{
		if ((mxGetM(prhs[n]) != N) || (mxGetN(prhs[n]) != N))
			mexErrMsgTxt("Dear student, the C_mat and K_mat matrices must be square with size identical to M_mat.");    //Improve this
	}

	if (mxGetN(prhs[3]) != 1)
		mexErrMsgTxt("Dear student, w_col must be a column vector.");
	size_t n_w_points=mxGetM(prhs[3]);

	if (mxGetM(prhs[4]) != 1)
		mexErrMsgTxt("Dear student, n_row must be a row vector.");
	size_t n_n_row=mxGetN(prhs[4]);

	if (mxGetM(prhs[5]) != 1)
		mexErrMsgTxt("Dear student, m_row must be a row vector.");

	if (mxGetN(prhs[5]) !=n_n_row)
		mexErrMsgTxt("Dear student, m_row & m_row must have the same size.");

	//Inputs
	Matrix<double> M_mat=ColMajor2RowMajor(N,N,mxGetPr(prhs[0]));
	Matrix<double> C_mat=ColMajor2RowMajor(N,N,mxGetPr(prhs[1]));
	Matrix<double> K_mat=ColMajor2RowMajor(N,N,mxGetPr(prhs[2]));
	Matrix<double> w_col=ColMajor2RowMajor(n_w_points,1,mxGetPr(prhs[3]));
	double* n_row=mxGetPr(prhs[4]);
	double* m_row=mxGetPr(prhs[5]);
	valarray<size_t> index_vec(n_n_row);
	for (unsigned int n=0; n<n_n_row; n++)
		index_vec[n]=(n_row[n]-1)*N+(m_row[n]-1); //This is the index_vec of n_row[n] row and m_row[n] column of a C_mat matrix (row major) and subscript starts with 0

	//Do Actual Computations
	Matrix<complex<double> > H_cols(n_w_points,n_n_row);
	try
	{
		H_cols=MDOF_FRF_Visc_Slow(M_mat,C_mat,K_mat,w_col,index_vec);
	}
	catch (char str[])
	{
		mexErrMsgTxt(str);
	}

	plhs[0] = mxCreateDoubleMatrix(n_w_points,n_n_row, mxCOMPLEX); // Create the output matrix
	double* Hn1_out_real=mxGetPr(plhs[0]);
	double* Hn1_out_imaginary=mxGetPi(plhs[0]);
	ComplexRowMajor2ColMajor(H_cols,Hn1_out_real,Hn1_out_imaginary);

	return;
}