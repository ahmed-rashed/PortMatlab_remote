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

Mat<complex<double> > MDOF_FRF_Visc(const Col<complex<double> >& EigValues_col, const Mat<complex<double> >& EigVectors_Normalized, const Mat<double>& w_column, const Row<double>& n_C, const Row<double>& m_C);

class MexFunction : public matlab::mex::Function
{
	ArrayFactory factory;
	const std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

public:
	void operator()(ArgumentList outputs, ArgumentList inputs)
	{
		size_t N, N_cols, N_w, N_col;
		checkArguments(outputs, inputs, N, N_cols, N_w, N_col);

		Col<complex<double> > EigValues_col(MatlabCol2Armadillo<complex<double> >(inputs[0]));
		Mat<complex<double> > EigVectors_Normalized(MatlabMat2Armadillo<complex<double> >(inputs[1]));
		Mat<double> w_column(MatlabMat2Armadillo<double>(inputs[2]));
		Row<double> n_C(MatlabRow2Armadillo<double>(inputs[3]));
		Row<double> m_C(MatlabRow2Armadillo<double>(inputs[4]));

		n_C=n_C-1;	//Convert from Matlab 1 based indexing to C 0 based indexing
		m_C=m_C-1;	//Convert from Matlab 1 based indexing to C 0 based indexing

		//Do Actual Computations
		Mat<complex<double> > H_w_n_m_cols(N_w,N_col);
		try
		{
			H_w_n_m_cols=MDOF_FRF_Visc(EigValues_col, EigVectors_Normalized, w_column, n_C, m_C);
		}
		catch (std::logic_error exceptt)
		{
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar(exceptt.what()) }));
		}

		//Outputs
		outputs[0] = factory.createArray<complex<double>>({ N_w, N_col }, H_w_n_m_cols.begin(), H_w_n_m_cols.end());	// Create the 1st output

		return;
	}

	void checkArguments(ArgumentList outputs, ArgumentList inputs, size_t &N, size_t &N_cols, size_t &N_w, size_t &N_col)
	{
		if (inputs.size() != 5)
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, This function needs 5 inputs.") }));

		if (outputs.size() > 1)
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, this function returns only one output.") }));

		N = inputs[1].getDimensions()[0];
		N_cols = inputs[1].getDimensions()[1];
		if (2 * N != N_cols)
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, The EigVectors_Normalized matrix must have columns twice as many as be rows.") }));

		if (inputs[0].getNumberOfElements() != 2 * N)
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, the EigValues_col vector must have length identical to the number of columns of EigVectors_Normalized matrix.") }));

		if (inputs[2].getDimensions()[1] != 1)
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, w_column must be a column vector.") }));
		N_w = inputs[2].getDimensions()[0];

		if (inputs[3].getDimensions()[0] != 1)
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, n_row must be a row vector.") }));
		N_col = inputs[3].getDimensions()[1];

		if ((inputs[4].getDimensions()[0] != 1) || (inputs[4].getDimensions()[1] != N_col))
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, m_row must have size identical to n_row.") }));
	};
};

