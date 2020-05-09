#include <complex>
#include <valarray>
#include "valarray_ComplexRealOperators.hpp"
#include "MexOperations.hpp"
#include "a_matrices.hpp"
#include "mex.hpp"
#include "mexAdapter.hpp"

using namespace std;
using namespace a_matrices;

using namespace matlab::data;
using matlab::mex::ArgumentList;

void MDOF_Eig_Visc(const Matrix<double>& M_mat,const Matrix<double>& C_mat,const Matrix<double>& K_mat, const bool& isPropotional, Matrix<complex<double> >& EigVectors_Normalized, valarray<complex<double> >& EigValues_vec);

class MexFunction : public matlab::mex::Function
{
	ArrayFactory factory;
	const std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

public:
	void operator()(ArgumentList outputs, ArgumentList inputs)
	{
		size_t N;
		checkArguments(outputs, inputs, N);

		Matrix<double> M_mat=ColMajor2RowMajor(static_cast<TypedArray<double>>(inputs[0]));
		Matrix<double> C_mat=ColMajor2RowMajor(static_cast<TypedArray<double>>(inputs[1]));
		Matrix<double> K_mat=ColMajor2RowMajor(static_cast<TypedArray<double>>(inputs[2]));

		bool isPropotional=false;
		if (inputs.size() == 4)
			isPropotional= inputs[3][0];

		//Do Actual Computations
		Matrix<complex<double> > EigVectors_Normalized(N,2*N);
		valarray<complex<double> > EigVal_vec(2*N);
		try
		{
			MDOF_Eig_Visc(M_mat,C_mat,K_mat,isPropotional, EigVectors_Normalized, EigVal_vec);
		}
		catch (char str[])
		{
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar(str) }));
		}

		//Outputs
		outputs[0] = std::move(RowMajor2ColMajor(EigVectors_Normalized));
		if (outputs.size() > 1)
			outputs[1] = factory.createArray<complex<double> >({1, 2 * N }, &EigVal_vec[0], &EigVal_vec[2*N]);

		return;
	}

	void checkArguments(ArgumentList outputs, ArgumentList inputs, size_t &N)
	{
		if ((inputs.size() != 3) && (inputs.size() != 4))
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, This function takes 3 or 4 inputs.") }));

		if ((outputs.size() != 1) && (outputs.size() != 2))
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, this function returns one or two outputs only.") }));

		N = inputs[0].getDimensions()[0];
		if (N != inputs[0].getDimensions()[1])
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, The M matrix must be square.") }));

		for (size_t n = 1; n < 3; n++)
		{
			if ((inputs[n].getDimensions()[0] != N) || (inputs[n].getDimensions()[1] != N))
				matlabPtr->feval(u"error", 0,
					std::vector<Array>({ factory.createScalar("Dear student, the C and K matrices must be square with size identical to M.") }));    //Improve this
		}

		if (inputs.size() == 4)
			if (inputs[3].getType() != ArrayType::LOGICAL || inputs[3].getNumberOfElements() != 1)
				matlabPtr->feval(u"error", 0,
					std::vector<Array>({ factory.createScalar("Dear student, isPropotional must be a logical scalar/integer.") }));
	};
};
