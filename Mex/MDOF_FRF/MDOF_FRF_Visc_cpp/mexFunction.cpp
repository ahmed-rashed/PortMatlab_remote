#include <complex>
#include <valarray>
#include "mex.h"
#include "valarray_ComplexRealOperators.h"
#include "MatrixOperations.h"
#include "a_matrices.h"

using namespace std;
using namespace a_matrices;

Matrix<complex<double> > MDOF_FRF_Visc(const valarray<complex<double> >& EigValues_vec, const Matrix<complex<double> >& EigVectors_Normalized, const Matrix<double>& w_column, const valarray<size_t>& n_vec, const valarray<size_t>& m_vec);

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 
    // Check for proper number of arguments
    if (nrhs != 5)
		mexErrMsgTxt("Dear student, This function needs 5 inputs."); 

	if (nlhs > 1)
		mexErrMsgTxt("Dear student, this function returns only one output."); 
    
	size_t N=mxGetM(prhs[1]);
	if (2*N != mxGetN(prhs[1]))
		mexErrMsgTxt("Dear student, The EigVectors_Normalized matrix must have columns twice as many as be rows.");

	if (mxGetNumberOfElements(prhs[0]) != 2*N)
		mexErrMsgTxt("Dear student, the EigValues_vec vector must have length identical to the number of columns of EigVectors_Normalized matrix.");

	if (mxGetN(prhs[2]) != 1)
		mexErrMsgTxt("Dear student, w_column must be a column vector.");
	size_t n_w=mxGetM(prhs[2]);

	if (mxGetM(prhs[3]) != 1)
		mexErrMsgTxt("Dear student, n_row must be a row vector.");
	size_t n_ii_row=mxGetN(prhs[3]);

	if ((mxGetM(prhs[4]) != 1) || (mxGetN(prhs[4]) != n_ii_row))
		mexErrMsgTxt("Dear student, m_row must have size identical to n_row.");

	//Inputs
	valarray<complex<double> > EigValues_vec=ComplexColMajor2RowMajor(1,2*N,mxGetPr(prhs[0]),mxGetPi(prhs[0])).va_copy();
	Matrix<complex<double> > EigVectors_Normalized=ComplexColMajor2RowMajor(N,2*N,mxGetPr(prhs[1]),mxGetPi(prhs[1]));
	Matrix<double> w_column=ColMajor2RowMajor(n_w,1,mxGetPr(prhs[2]));
	double* ii_matlab=mxGetPr(prhs[3]);
	double* jj_matlab=mxGetPr(prhs[4]);

	valarray<size_t> n_vec(n_ii_row);
	valarray<size_t> m_vec(n_ii_row);
	for (unsigned int n=0; n<n_ii_row; n++)
	{
		n_vec[n]=ii_matlab[n]-1;	//Convert from Matlab 1 based indexing to C 0 based indexing
		m_vec[n]=jj_matlab[n]-1;	//Convert from Matlab 1 based indexing to C 0 based indexing
	}

	//Do Actual Computations
	Matrix<complex<double> > H_w_n_m_cols(n_w,n_ii_row);
	try
	{
		H_w_n_m_cols=MDOF_FRF_Visc(EigValues_vec, EigVectors_Normalized, w_column, n_vec, m_vec);
	}
	catch (char str[])
	{
		mexErrMsgTxt(str);
	}

	plhs[0] = mxCreateDoubleMatrix(n_w,n_ii_row, mxCOMPLEX); // Create the output matrix
	double* H_w_mat_out_real=mxGetPr(plhs[0]);
	double* H_w_mat_out_imaginary=mxGetPi(plhs[0]);
	ComplexRowMajor2ColMajor(H_w_n_m_cols,H_w_mat_out_real,H_w_mat_out_imaginary);

	return;
}
