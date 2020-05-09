#include <valarray>
#include "mex.hpp"
#include "mexAdapter.hpp"

using namespace matlab::data;
using matlab::mex::ArgumentList;

std::valarray<double> Free_Response(const double& w_n, const double& zeta, const double& x0, const double& v0, const std::valarray<double>& t);

class MexFunction : public matlab::mex::Function
{
	ArrayFactory factory;
	
public:
	void operator()(ArgumentList outputs, ArgumentList inputs)
	{
		checkArguments(outputs, inputs);
		
		const double w_n=inputs[0][0];
		const double zeta= inputs[1][0];
		const double x0= inputs[2][0];
		const double v0= inputs[3][0];
		const size_t N= inputs[4].getNumberOfElements();
		std::valarray<double> t(N);
		for (size_t n = 0; n < N; n++)
			t[n]= inputs[4][n];
	
		// Do the actual computations
		std::valarray<double> x=Free_Response(w_n,zeta,x0,v0,t);

		TypedArray<double> x_out=
		outputs[0] = factory.createArray<double>(inputs[4].getDimensions(), &x[0], &x[N]);

		return;
	}

	void checkArguments(ArgumentList outputs, ArgumentList inputs)
	{
		const std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

		if (inputs.size() != 5)
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, This function needs 5 inputs.") }));

		if (outputs.size() > 1)
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, this function returns only one output.") }));

		for (size_t n = 0; n < 4; n++)
		{
			if (inputs[n].getNumberOfElements() != 1)
				matlabPtr->feval(u"error", 0,
					std::vector<Array>({ factory.createScalar("Dear student, first four inputs must be scalars.") }));    //Improve this
		}
	};
};