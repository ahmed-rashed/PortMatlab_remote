#ifndef MexOperations_hpp
#define MexOperations_hpp

#include "a_matrices.hpp"
#include "MatlabDataArray/TypedArrayRef.hpp"

template <class T>
a_matrices::Matrix<T> MatlabMat2Matrix(const matlab::data::Array& matlabMatrixColMajor)
{
	size_t m = matlabMatrixColMajor.getDimensions()[0];
	size_t n = matlabMatrixColMajor.getDimensions()[1];

	const TypedArray<T> matlabTypedArrayColMajor = matlabMatrixColMajor;
	a_matrices::Matrix<T> ResRowMajor(n, m);
	size_t ii = 0;
	for (auto& elem : matlabTypedArrayColMajor)
	{
		ResRowMajor(ii) = elem;
		ii++;
	}

	return transpose(ResRowMajor);
};

template <class T>
std::valarray<T> MatlabVec2valarray(const matlab::data::Array& matlabVector)
{
	size_t N = matlabVector.getNumberOfElements();

	const TypedArray<T> matlabTypedArray = matlabVector;
	std::valarray<T> Res(N);
	size_t ii = 0;
	for (auto& elem : matlabTypedArray)
	{
		Res[ii] = elem;
		ii++;
	}

	return Res;
};

template <class T>
matlab::data::TypedArray<T> RowMajor2ColMajor(const a_matrices::Matrix<T>& mat_RowMajor)
{
	a_matrices::Matrix<T> mat_temp=transpose(mat_RowMajor);

	size_t m = mat_RowMajor.rows;
	size_t n = mat_RowMajor.cols;
	matlab::data::ArrayFactory factory;
	matlab::data::TypedArray<T> matlabMatrixColMajor = factory.createArray<T>({ m, n }, mat_temp.begin(), mat_temp.end());

	return matlabMatrixColMajor;
};

#endif