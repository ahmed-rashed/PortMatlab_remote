#include "mex.h"
void Free_Response(const double* const p_w_n, const double* const p_zeta, const double* const p_x0, const double* const p_v0, const double t[], double x[], const size_t* const p_N);
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 
    unsigned int n;
	double w_n;
	double zeta;
	double x0;
	double v0;
	size_t N;
	double *t; 
	double *x;

	// Check for proper number of arguments
    if (nrhs != 5)
		mexErrMsgTxt("Dear student, This function needs 5 inputs."); 

	if (nlhs > 1)
		mexErrMsgTxt("Dear student, this function returns only one output."); 
    
	for (n=0; n<4; n++)
	{
		if (mxGetNumberOfElements(prhs[n]) != 1)
			mexErrMsgTxt("Dear student, first four inputs must be scalars.");    //Improve this
	}

	w_n=mxGetScalar(prhs[0]);
	zeta=mxGetScalar(prhs[1]);
	x0=mxGetScalar(prhs[2]);
	v0=mxGetScalar(prhs[3]);
	N=mxGetNumberOfElements(prhs[4]);
	t = mxGetPr(prhs[4]); 
	plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[4]),mxGetN(prhs[4]), mxREAL); // Create the output matrix with the same size as t
	x = mxGetPr(plhs[0]);

	// Do the actual computations
	Free_Response(&w_n , &zeta , &x0 , &v0 , t , x , &N);

	return;
}
