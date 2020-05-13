#include <complex>

#include "ArmadilloMatrixOperations.hpp"
#include "ArmadilloMexOperations.hpp"

#include "mex.hpp"
#include "mexAdapter.hpp"	//Include mexAdapter.hpp only once with the MexFunction class definition in MEX applications that span multiple files.
							//[www.mathworks.com/help/matlab/matlab_external/structure-of-c-mex-function.html]

using namespace matlab::data;
using matlab::mex::ArgumentList;

using namespace std;
using namespace arma;

void MDOF_Eig_Visc(const Mat<double>& M_mat,const Mat<double>& C_mat,const Mat<double>& K_mat, const bool& isPropotional, Mat<complex<double> >& EigVectors_Normalized, Col<complex<double> >& EigValues_col);

class MexFunction : public matlab::mex::Function
{
	ArrayFactory factory;
	const std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

public:
	void operator()(ArgumentList outputs, ArgumentList inputs)
	{
		size_t N;
		checkArguments(outputs, inputs, N);

		Mat<double> M_mat(MatlabMat2Armadillo<double>(inputs[0]));
		Mat<double> C_mat(MatlabMat2Armadillo<double>(inputs[1]));
		Mat<double> K_mat(MatlabMat2Armadillo<double>(inputs[2]));

		bool isPropotional=false;
		if (inputs.size() == 4)
			isPropotional= inputs[3][0];

		//Do Actual Computations
		Mat<complex<double> > EigVectors_Normalized(N,2*N);
		Col<complex<double> > EigValues_col(2*N);
		try
		{
			MDOF_Eig_Visc(M_mat,C_mat,K_mat,isPropotional, EigVectors_Normalized, EigValues_col);
		}
		catch (std::logic_error exceptt)
		{
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar(exceptt.what()) }));
		}

		//Outputs
		outputs[0] = factory.createArray<complex<double>>({ N,2 * N }, EigVectors_Normalized.begin(), EigVectors_Normalized.end());	// Create the 1st output
		if (outputs.size() > 1)	// Create the 2nd output
			outputs[1] = factory.createArray<complex<double>>({ 2 * N,1 }, EigValues_col.begin(), EigValues_col.end());

		return;
	}

	void checkArguments(ArgumentList outputs, ArgumentList inputs, size_t& N)
	{
		if ((inputs.size() != 3) && (inputs.size() != 4))
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, This function takes 3 or 4 inputs.") }));

		if ((outputs.size() != 1) && (outputs.size() != 2))
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, this function returns one or two outputs only.") }));

		N = inputs[0].getDimensions()[0];
		size_t n_cols = inputs[0].getDimensions()[1];
		if (N != n_cols)
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

