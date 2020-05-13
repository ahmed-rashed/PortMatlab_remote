#include <complex>
#include <valarray>
#include <math.h>
#include "a_matrices.hpp"
#include "MatrixOperations.hpp"
#include "MexOperations.hpp"
#include "mex.hpp"
#include "mexAdapter.hpp"	//Include mexAdapter.hpp only once with the MexFunction class definition in MEX applications that span multiple files.
							//[www.mathworks.com/help/matlab/matlab_external/structure-of-c-mex-function.html]

using namespace std;
using namespace a_matrices;

using namespace matlab::data;
using matlab::mex::ArgumentList;

Matrix<complex<double> > MDOF_FRF_Visc_Slow(const Matrix<double>& M_mat,const Matrix<double>& C_mat,const Matrix<double>& K_mat,const Matrix<double>& w_col, const valarray<size_t>& index_vec);

class MexFunction : public matlab::mex::Function
{
	ArrayFactory factory;

public:
	void operator()(ArgumentList outputs, ArgumentList inputs)
	{
		size_t N, N_w_points, N_n_row;
		checkArguments(outputs, inputs, N, N_w_points, N_n_row);

		Matrix<double> M_mat(MatlabMat2Matrix<double>(inputs[0]));
		Matrix<double> C_mat(MatlabMat2Matrix<double>(inputs[1]));
		Matrix<double> K_mat(MatlabMat2Matrix<double>(inputs[2]));
		Matrix<double> w_col(MatlabMat2Matrix<double>(inputs[3]));
		valarray<double> ii_vec_temp(MatlabVec2valarray<double>(inputs[4]));
		valarray<double> jj_vec_temp(MatlabVec2valarray<double>(inputs[5]));

		valarray<size_t> ii_vec(ii_vec_temp.size()), jj_vec(jj_vec_temp.size());
		copy(begin(ii_vec_temp), end(ii_vec_temp), begin(ii_vec));
		copy(begin(jj_vec_temp), end(jj_vec_temp), begin(jj_vec));
		valarray<size_t> index_vec = RowMajorSub2Ind(N, ii_vec-1, jj_vec-1); //Note that ii_row and jj_row and 1-based Matlab subscript

		//Do Actual Computations
		Matrix<complex<double> > H_cols(N_w_points,N_n_row);
		try
		{
			H_cols=MDOF_FRF_Visc_Slow(M_mat,C_mat,K_mat,w_col,index_vec);
		}
		catch (char str[])
		{
			const std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar(str) }));
		}

		outputs[0] = std::move(RowMajor2ColMajor(H_cols));

		return;
	}

	void checkArguments(ArgumentList outputs, ArgumentList inputs, size_t &N, size_t &N_w_points, size_t &N_n_row)
	{
		const std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

		if (inputs.size() != 6)
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, This function needs 6 inputs.") }));

		if (outputs.size() > 1)
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, this function returns only one output.") }));

		N = inputs[0].getDimensions()[0];
		if (inputs[0].getType() != ArrayType::DOUBLE ||
			inputs[0].getDimensions()[1] != N)
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, The M_mat matrix must be square.") }));

		for (size_t n = 1; n < 3; n++)
			if ((inputs[n].getType() != ArrayType::DOUBLE) ||
				(inputs[n].getDimensions()[0] != N) ||
				(inputs[n].getDimensions()[1] != N))
				matlabPtr->feval(u"error", 0,
					std::vector<Array>({ factory.createScalar("Dear student, the C_mat and K_mat matrices must be square with size identical to M_mat.") }));    //Improve this

		if (inputs[3].getType() != ArrayType::DOUBLE || 
			inputs[3].getDimensions()[1] != 1)
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, w_col must be a column vector.") }));
		N_w_points = inputs[3].getDimensions()[0];

		if (inputs[4].getType() != ArrayType::DOUBLE || 
			inputs[4].getDimensions()[0] != 1)
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, n_row must be a row vector.") }));
		N_n_row = inputs[4].getDimensions()[1];

		if (inputs[5].getType() != ArrayType::DOUBLE || 
			inputs[5].getDimensions()[0] != 1)
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, m_row must be a row vector.") }));

		if (inputs[5].getDimensions()[1] != N_n_row)
			matlabPtr->feval(u"error", 0,
				std::vector<Array>({ factory.createScalar("Dear student, m_row & m_row must have the same size.") }));
	};
};