#include <complex>
#include <valarray>
#include <math.h>
#include "a_matrices.hpp"
#include "MexOperations.hpp"
#include "mex.hpp"
#include "mexAdapter.hpp"

using namespace std;
using namespace a_matrices;

using namespace matlab::data;
using matlab::mex::ArgumentList;

Matrix<complex<double> > MDOF_FRF_Visc_Slow(const Matrix<double>& M_mat,const Matrix<double>& C_mat,const Matrix<double>& K_mat,const Matrix<double>& w_col, const valarray<size_t>& index_vec);

class MexFunction : public matlab::mex::Function
{
	ArrayFactory factory;
	const std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

public:
	void operator()(ArgumentList outputs, ArgumentList inputs)
	{
		size_t N, N_w_points, N_n_row;
		checkArguments(outputs, inputs, N, N_w_points, N_n_row);

		/*Matrix<double> M_mat=ColMajor2RowMajor(static_cast<TypedArray<double> >(inputs[0]));
		Matrix<double> C_mat=ColMajor2RowMajor(static_cast<TypedArray<double> >(inputs[1]));
		Matrix<double> K_mat=ColMajor2RowMajor(static_cast<TypedArray<double> >(inputs[2]));
		Matrix<double> w_col=ColMajor2RowMajor(static_cast<TypedArray<double> >(inputs[3]));

		valarray<size_t> ii_vec(N_n_row), jj_vec(N_n_row);
		for (size_t n = 0; n < N_n_row; n++)
		{
			ii_vec[n] = static_cast<size_t>(inputs[4][n]);
			jj_vec[n] = static_cast<size_t>(inputs[5][n]);
		}*/

		//valarray<size_t> index_vec= (ii_vec - 1) * N + (jj_vec - 1); //index_vec is a row-major and 0-based C index, while ii_row and jj_row and 1-based Matlab subscript

		//Do Actual Computations
		Matrix<complex<double> > H_cols(N_w_points,N_n_row);
		try
		{
			//H_cols=MDOF_FRF_Visc_Slow(M_mat,C_mat,K_mat,w_col,index_vec);
		}
		catch (char str[])
		{
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar(str) }));
		}

		outputs[0] = std::move(RowMajor2ColMajor(H_cols));

		return;
	}

	void checkArguments(ArgumentList outputs, ArgumentList inputs, size_t &N, size_t &N_w_points, size_t &N_n_row)
	{
		if (inputs.size() != 6)
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, This function needs 6 inputs.") }));

		if (outputs.size() > 1)
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, this function returns only one output.") }));

		N = inputs[0].getDimensions()[0];
		if (N != inputs[0].getDimensions()[1])
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, The M_mat matrix must be square.") }));

		for (size_t n = 1; n < 3; n++)
		{
			if ((inputs[n].getDimensions()[0] != N) || (inputs[n].getDimensions()[1] != N))
				matlabPtr->feval(u"error", 0,
					std::vector<Array>({ factory.createScalar("Dear student, the C_mat and K_mat matrices must be square with size identical to M_mat.") }));    //Improve this
		}

		if (inputs[3].getDimensions()[1] != 1)
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, w_col must be a column vector.") }));
		N_w_points = inputs[3].getDimensions()[0];

		if (inputs[4].getDimensions()[0] != 1)
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, n_row must be a row vector.") }));
		N_n_row = inputs[4].getDimensions()[1];

		if (inputs[5].getDimensions()[0] != 1)
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, m_row must be a row vector.") }));

		if (inputs[5].getDimensions()[1] != N_n_row)
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, m_row & m_row must have the same size.") }));
	};
};