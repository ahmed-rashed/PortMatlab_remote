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

Matrix<complex<double> > MDOF_FRF_Visc(const valarray<complex<double> >& EigValues_vec, const Matrix<complex<double> >& EigVectors_Normalized, const Matrix<double>& w_column, const valarray<size_t>& n_vec, const valarray<size_t>& m_vec);

class MexFunction : public matlab::mex::Function
{
	ArrayFactory factory;
	const std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

public:
	void operator()(ArgumentList outputs, ArgumentList inputs)
	{
		size_t N, N_ii_row;
		checkArguments(outputs, inputs, N, N_ii_row);

		valarray<complex<double> > EigValues_vec(N);
		for (size_t n = 0; n < N; n++)
			EigValues_vec[n] = inputs[0][n];

		Matrix<complex<double> > EigVectors_Normalized = ColMajor2RowMajor(static_cast<TypedArray<complex<double> >>(inputs[1]));

		size_t N_w = inputs[2].getDimensions()[0];
		Matrix<double> w_column= ColMajor2RowMajor(static_cast<TypedArray<double> >(inputs[2]));

		valarray<size_t> n_vec(N_ii_row), m_vec(N_ii_row);
		for (size_t n = 0; n < N_ii_row; n++)
		{
			n_vec[n] = inputs[3][n];
			m_vec[n] = inputs[4][n];
		}
		n_vec = n_vec - 1;	//Convert from Matlab 1-based indexing to C 0-based indexing
		m_vec = m_vec - 1;	//Convert from Matlab 1-based indexing to C 0-based indexing

		//Do Actual Computations
		Matrix<complex<double> > H_w_n_m_cols(N_w, N_ii_row);
		try
		{
			H_w_n_m_cols = MDOF_FRF_Visc(EigValues_vec, EigVectors_Normalized, w_column, n_vec, m_vec);
		}
		catch (char str[])
		{
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar(str) }));
		}

		outputs[0] = RowMajor2ColMajor(H_w_n_m_cols);

		return;
	}

	void checkArguments(ArgumentList outputs, ArgumentList inputs, size_t& N, size_t& N_ii_row)
	{
		if (inputs.size() != 5)
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, This function needs 5 inputs.") }));

		if (outputs.size() > 1)
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, this function returns only one output.") }));

		N = inputs[1].getDimensions()[0];
		if (2 * N != inputs[1].getDimensions()[1])
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, The EigVectors_Normalized matrix must have columns twice as many as be rows.") }));

		if (inputs[0].getNumberOfElements() != 2 * N)
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, the EigValues_vec vector must have length identical to the number of columns of EigVectors_Normalized matrix.") }));

		if (inputs[2].getDimensions()[1] != 1)
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, w_column must be a column vector.") }));

		if (inputs[3].getDimensions()[0] != 1)
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, n_row must be a row vector.") }));

		N_ii_row = inputs[3].getDimensions()[1];
		if ((inputs[4].getDimensions()[0] != 1) || (inputs[4].getDimensions()[1] != N_ii_row))
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, m_row must have size identical to n_row.") }));
	};
};
