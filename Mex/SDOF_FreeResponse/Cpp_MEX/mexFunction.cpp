#include <cstring>
#include <valarray>
#include "mex.h"

std::valarray<double> Free_Response(const double& w_n, const double& zeta, const double& x0, const double& v0, const std::valarray<double>& t);

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 
    // Check for proper number of arguments
    if (nrhs != 5)
		mexErrMsgTxt("Dear student, This function needs 5 inputs."); 

	if (nlhs > 1)
		mexErrMsgTxt("Dear student, this function returns only one output."); 
    
	for (unsigned int n=0; n<4; n++)
	{
		if (mxGetNumberOfElements(prhs[n]) != 1)
			mexErrMsgTxt("Dear student, first four inputs must be scalars.");    //Improve this
	}

	double w_n=mxGetScalar(prhs[0]);
	double zeta=mxGetScalar(prhs[1]);
	double x0=mxGetScalar(prhs[2]);
	double v0=mxGetScalar(prhs[3]);
	size_t N=mxGetNumberOfElements(prhs[4]);
	//double *t = mxGetPr(prhs[4]);
	std::valarray<double> t(mxGetPr(prhs[4]),N);
	
	// Do the actual computations
	std::valarray<double> x=Free_Response(w_n,zeta,x0,v0,t);

	plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[4]),mxGetN(prhs[4]), mxREAL); // Create the output matrix with the same size as t
	double* x_out=mxGetPr(plhs[0]);
	for (unsigned int n=0; n<N; n++)
	{
		x_out[n]=x[n];
	}
	//memcpy(mxGetPr(plhs[0]), &x[0], N*mxGetElementSize(plhs[0]));
	memcpy(mxGetPr(plhs[0]), &x[0], N*sizeof(x[0]));

	return;
}
