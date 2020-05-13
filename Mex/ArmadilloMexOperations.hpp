#ifndef ArmadilloMexOperations_hpp
#define ArmadilloMexOperations_hpp

#include <armadillo>
#include "MatlabDataArray/TypedArrayRef.hpp"

using namespace std;
using namespace arma;

template <class T>
Mat<T> MatlabMat2Armadillo(const matlab::data::Array& matlabMatrix)
{
	size_t m = matlabMatrix.getDimensions()[0];
	size_t n = matlabMatrix.getDimensions()[1];

	const TypedArray<T> matlabTypedArray = matlabMatrix;
	Mat<T> ResMat(m, n);
	size_t ii = 0;
	for (auto& elem : matlabTypedArray)
	{
		ResMat(ii) = elem;
		ii++;
	}

	return ResMat;
};

template <class T>
Col<T> MatlabCol2Armadillo(const matlab::data::Array& matlabCol)
{
	size_t m = matlabCol.getDimensions()[0];

	const TypedArray<T> matlabTypedArray = matlabCol;
	Col<T> ResCol(m);
	size_t ii = 0;
	for (auto& elem : matlabTypedArray)
	{
		ResCol(ii) = elem;
		ii++;
	}

	return ResCol;
};

template <class T>
Row<T> MatlabRow2Armadillo(const matlab::data::Array& matlabRow)
{
	size_t n = matlabRow.getDimensions()[1];

	const TypedArray<T> matlabTypedArray = matlabRow;
	Row<T> ResRow(n);
	size_t ii = 0;
	for (auto& elem : matlabTypedArray)
	{
		ResRow(ii) = elem;
		ii++;
	}

	return ResRow;
};

#endif