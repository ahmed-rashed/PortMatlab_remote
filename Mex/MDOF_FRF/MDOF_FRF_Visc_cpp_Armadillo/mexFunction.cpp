#include <complex>

#include "mex.h"

#include "ArmadilloMatrixOperations.h"
#include "armaMex.hpp"

using namespace std;
using namespace arma;

Mat<complex<double> > MDOF_FRF_Visc(const Col<complex<double> >& EigValues_col, const Mat<complex<double> >& EigVectors_Normalized, const Mat<double>& w_column, const Row<double>& n_C, const Row<double>& m_C);

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 
    // Check for proper number of arguments
    if (nrhs != 5)
		mexErrMsgTxt("Dear student, This function needs 5 inputs."); 

	if (nlhs > 1)
		mexErrMsgTxt("Dear student, this function returns only one output."); 
    
	size_t N=mxGetM(prhs[1]);
	size_t n_cols=mxGetN(prhs[1]);
	if (2*N != n_cols)
		mexErrMsgTxt("Dear student, The EigVectors_Normalized matrix must have columns twice as many as be rows.");

	if (mxGetNumberOfElements(prhs[0]) != 2*N)
		mexErrMsgTxt("Dear student, the EigValues_col vector must have length identical to the number of columns of EigVectors_Normalized matrix.");

	if (mxGetN(prhs[2]) != 1)
		mexErrMsgTxt("Dear student, w_column must be a column vector.");
	size_t n_w=mxGetM(prhs[2]);

	if (mxGetM(prhs[3]) != 1)
		mexErrMsgTxt("Dear student, n_row must be a row vector.");
	size_t n_col=mxGetN(prhs[3]);

	if ((mxGetM(prhs[4]) != 1) || (mxGetN(prhs[4]) != n_col))
		mexErrMsgTxt("Dear student, m_row must have size identical to n_row.");

	//Inputs
	Col<complex<double> > EigValues_col=armaGetCx(prhs[0]);
	Mat<complex<double> > EigVectors_Normalized=armaGetCx(prhs[1]);
	Mat<double> w_column=armaGetPr(prhs[2]);
	Row<double> n_C=armaGetPr(prhs[3]);
	Row<double> m_C=armaGetPr(prhs[4]);

	n_C=n_C-1;	//Convert from Matlab 1 based indexing to C 0 based indexing
	m_C=m_C-1;	//Convert from Matlab 1 based indexing to C 0 based indexing

	//Do Actual Computations
	Mat<complex<double> > H_w_n_m_cols(n_w,n_col);
	try
	{
		H_w_n_m_cols=MDOF_FRF_Visc(EigValues_col, EigVectors_Normalized, w_column, n_C, m_C);
	}
	catch (std::logic_error exceptt)
	{
		mexErrMsgTxt(exceptt.what());
	}

	//Outputs
	plhs[0] = armaCreateMxMatrix(n_w,n_col,mxDOUBLE_CLASS,mxCOMPLEX);	// Create the output matrix
	armaSetCx(plhs[0], H_w_n_m_cols);

	return;
}
