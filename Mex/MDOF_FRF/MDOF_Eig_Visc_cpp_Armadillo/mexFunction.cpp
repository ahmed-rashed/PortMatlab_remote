#include <complex>
#include "mex.h"

#include "ArmadilloMatrixOperations.h"
#include "armaMex.hpp"

using namespace std;
using namespace arma;

void MDOF_Eig_Visc(const Mat<double>& M_mat,const Mat<double>& C_mat,const Mat<double>& K_mat, const bool& isPropotional, Mat<complex<double> >& EigVectors_Normalized, Col<complex<double> >& EigValues_col);

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 
    // Check for proper number of arguments
    if ((nrhs != 3) && (nrhs != 4))
		mexErrMsgTxt("Dear student, This function takes 3 or 4 inputs."); 

	if ((nlhs != 1) && (nlhs != 2))
		mexErrMsgTxt("Dear student, this function returns one or two outputs only.");
    
	size_t N=mxGetM(prhs[0]);
	size_t n_cols=mxGetN(prhs[0]);
	if (N != n_cols)
		mexErrMsgTxt("Dear student, The M matrix must be square.");

	for (unsigned int n=1; n<3; n++)
	{
		if ((mxGetM(prhs[n]) != N) || (mxGetN(prhs[n]) != N))
			mexErrMsgTxt("Dear student, the C and K matrices must be square with size identical to M.");    //Improve this
	}

	if (nrhs == 4)
		if (!mxIsLogicalScalar(prhs[3]))
			mexErrMsgTxt("Dear student, isPropotional must be a logical scalar/integer.");

	//Inputs
	Mat<double> M_mat=armaGetPr(prhs[0]);
	Mat<double> C_mat=armaGetPr(prhs[1]);
	Mat<double> K_mat=armaGetPr(prhs[2]);

	bool isPropotional=false;
	if (nrhs == 4)
		isPropotional=armaGetScalar<bool>(prhs[3]);

	//Do Actual Computations
	Mat<complex<double> > EigVectors_Normalized(N,2*N);
	Col<complex<double> > EigValues_col(2*N);
	try
	{
		MDOF_Eig_Visc(M_mat,C_mat,K_mat,isPropotional, EigVectors_Normalized, EigValues_col);
	}
	catch (std::logic_error exceptt)
	{
		mexErrMsgTxt(exceptt.what());
	}

	//Outputs
	plhs[0] = armaCreateMxMatrix(N,2*N,mxDOUBLE_CLASS,mxCOMPLEX);	// Create the 1st output matrix
	armaSetCx(plhs[0], EigVectors_Normalized);

	if (nlhs > 1)
	{
		plhs[1] = armaCreateMxMatrix(2*N,1,mxDOUBLE_CLASS,mxCOMPLEX);	// Create the 2nd output matrix
		armaSetCx(plhs[1],EigValues_col);
	}

	return;
}
