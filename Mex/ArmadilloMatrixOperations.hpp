#ifndef ArmadilloMatrixOperations_hpp
#define ArmadilloMatrixOperations_hpp

#include <armadillo>

using namespace std;
using namespace arma;


inline Col<int> stepvec(int nStart, int step, int nEnd)
{
	size_t N=floor((static_cast<double>(nEnd)-nStart)/step)+1;
	Col<int> outVec(N);
	if (N>=1)
	{
		outVec[0]=nStart;
		for (int n=1; n<N; n++)
			outVec[n]=outVec[n-1]+step;
	}

	return outVec;
}

template <class T>
inline uword ColMajorSub2Ind_Scalar(const size_t& n_rows, const T& row, const T& col)
{
	return static_cast<uword>(col*n_rows+row);
}

//The following function is useless, since Armadillo is column major
template <class T>
inline uword RowMajorSub2Ind_Scalar(const size_t& n_cols, const T& row, const T& col)
{
	return ColMajorSub2Ind_Scalar(n_cols,  col, row);
}

template <class T>
inline void ColMajorInd2Sub_Scalar(const size_t& n_rows, const T& ind, uword& row, uword& col)
{
	ldiv_t temp=div(static_cast<long>(ind),n_rows);
	row=temp.rem;
	col=temp.quot;
}

//The following function is useless, since Armadillo is column major
template <class T>
inline void RowMajorInd2Sub_Scalar(const size_t& n_cols, const T& ind, uword& row, uword& col)
{
	ColMajorInd2Sub_Scalar(n_cols, ind, col, row);
}



template <class T>
inline Mat<uword> ColMajorSub2Ind(const size_t& n_rows, const Mat<T>& row, const Mat<T>& col)
{
	return conv_to<Mat<uword> >::from(col*n_rows+row);
}

//The following function is useless, since Armadillo is column major
template <class T>
inline Mat<uword> RowMajorSub2Ind(const size_t& n_cols, const Mat<T>& row, const Mat<T>& col)
{
	return ColMajorSub2Ind(n_cols,  col, row);
}

template <class T>
inline void ColMajorInd2Sub(const size_t& n_rows, const Mat<T>& ind, Mat<uword>& row, Mat<uword>& col)
{
	for (unsigned int n=0; n<ind.size(); n++)
	{
		ColMajorInd2Sub_Scalar(n_rows, ind[n], row[n], col[n]);
	}
}

//The following function is useless, since Armadillo is column major
template <class T>
inline void RowMajorInd2Sub(const size_t& n_cols, const Mat<T>& ind, Mat<uword>& row, Mat<uword>& col)
{
	ColMajorInd2Sub(n_cols, ind, col, row);
}


#endif