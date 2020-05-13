#ifndef MexOperations_hpp
#define MexOperations_hpp

#include "a_matrices.hpp"
#include "MatlabDataArray/TypedArrayRef.hpp"

template <class T>
a_matrices::Matrix<T> ColMajor2RowMajor(const matlab::data::TypedArray<T>& matlabArrayColMajor)
{
	size_t m = matlabArrayColMajor.getDimensions()[0];
	size_t n = matlabArrayColMajor.getDimensions()[1];
	a_matrices::Matrix<T> ResRowMajor(m, n);

	for (size_t i = 0; i < m; i++)
		for (size_t j = 0; j < n; j++)
			ResRowMajor[i][j] = matlabArrayColMajor[j * m + i];

	return ResRowMajor;
};

template <class T>
matlab::data::TypedArray<T> RowMajor2ColMajor(const a_matrices::Matrix<T>& H_RowMajor)
{
	size_t m = H_RowMajor.rows;
	size_t n = H_RowMajor.cols;
	matlab::data::ArrayFactory factory;
	matlab::data::TypedArray<T> matlabArrayColMajor = factory.createArray<T>({ m, n });

	for (size_t j = 0; j < n; j++)
		for (size_t i = 0; i < m; i++)
			matlabArrayColMajor[j * m + i] = H_RowMajor[i][j];

	return matlabArrayColMajor;
};

#endif