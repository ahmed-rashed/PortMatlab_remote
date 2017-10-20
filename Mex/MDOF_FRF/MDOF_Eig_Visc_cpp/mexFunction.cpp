#include <complex>
#include <valarray>
#include "mex.h"
#include "valarray_ComplexRealOperators.h"
#include "MatrixOperations.h"
#include "a_matrices.h"

using namespace std;
using namespace a_matrices;

void MDOF_Eig_Visc(const Matrix<double>& M_mat,const Matrix<double>& C_mat,const Matrix<double>& K_mat, const bool& isPropotional, Matrix<complex<double> >& EigVectors_Normalized, valarray<complex<double> >& EigValues_vec);

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 
    // Check for proper number of arguments
    if ((nrhs != 3) && (nrhs != 4))
		mexErrMsgTxt("Dear student, This function takes 3 or 4 inputs."); 

	if ((nlhs != 1) && (nlhs != 2))
		mexErrMsgTxt("Dear student, this function returns one or two outputs only.");
    
	size_t N=mxGetM(prhs[0]);
	if (N != mxGetN(prhs[0]))
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
	Matrix<double> M_mat=ColMajor2RowMajor(N,N,mxGetPr(prhs[0]));
	Matrix<double> C_mat=ColMajor2RowMajor(N,N,mxGetPr(prhs[1]));
	Matrix<double> K_mat=ColMajor2RowMajor(N,N,mxGetPr(prhs[2]));

	bool isPropotional=false;
	if (nrhs == 4)
		isPropotional=mxIsLogicalScalarTrue(prhs[3]);

	//Do Actual Computations
	Matrix<complex<double> > EigVectors_Normalized(N,2*N);
	valarray<complex<double> > EigVal_vec(2*N);
	try
	{
		MDOF_Eig_Visc(M_mat,C_mat,K_mat,isPropotional, EigVectors_Normalized, EigVal_vec);
	}
	catch (char str[])
	{
		mexErrMsgTxt(str);
	}

	//Outputs
	plhs[0] = mxCreateDoubleMatrix(N,2*N, mxCOMPLEX); // Create the 1st output matrix
	double* EigVectors_Normalized_out_real_p=mxGetPr(plhs[0]);
	double* EigVectors_Normalized_out_imaginary_p=mxGetPi(plhs[0]);
	ComplexRowMajor2ColMajor(EigVectors_Normalized,EigVectors_Normalized_out_real_p,EigVectors_Normalized_out_imaginary_p);

	if (nlhs > 1)
	{
		plhs[1] = mxCreateDoubleMatrix(2*N,1, mxCOMPLEX); // Create the 2nd output matrix
		double* EigValues_mat_out_real_p=mxGetPr(plhs[1]);
		double* EigValues_mat_out_imaginary_p=mxGetPi(plhs[1]);

		for (size_t i=0; i<2*N;i++)
		{
			EigValues_mat_out_real_p[i]=EigVal_vec[i].real();
			EigValues_mat_out_imaginary_p[i]=EigVal_vec[i].imag();
		}
	}

	return;
}
