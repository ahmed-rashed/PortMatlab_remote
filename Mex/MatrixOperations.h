#ifndef MatrixOperations_h
#define MatrixOperations_h

#include <valarray>
#include <complex>
#include "a_matrices.h"

using namespace std;
using namespace a_matrices;

inline size_t RowMajorSub2Ind(const size_t& n_cols, const size_t& row, const size_t& col)
{
	return row*n_cols+col;
}

inline size_t ColMajorSub2Ind(const size_t& n_rows, const size_t& row, const size_t& col)
{
	return RowMajorSub2Ind(n_rows, col, row);
}

inline void RowMajorInd2Sub(const size_t& n_cols, const size_t& ind, size_t& row, size_t& col)
{
	ldiv_t temp=div(static_cast<long>(ind),n_cols);
	row=temp.quot;
	col=temp.rem;
}

inline void ColMajorInd2Sub(const size_t& n_rows, const size_t& ind, size_t& row, size_t& col)
{
	RowMajorInd2Sub(n_rows, ind, col, row);
}


template <template<typename> class T>
inline T<size_t> RowMajorSub2Ind(const size_t& n_cols, const T<size_t>& row, const T<size_t>& col)
{
	return row*n_cols+col;
}

template <template<typename> class T>
inline T<size_t> ColMajorSub2Ind(const size_t& n_rows, const T<size_t>& row, const T<size_t>& col)
{
	return RowMajorSub2Ind(n_rows, col, row);
}

template <template<typename> class T>
inline void RowMajorInd2Sub(const size_t& n_cols, const T<size_t>& ind, T<size_t>& row, T<size_t>& col)
{
	for (unsigned int n=0; n<ind.size(); n++)
	{
		RowMajorInd2Sub(n_cols, ind[n], row[n], col[n]);
	}
}

template <template<typename> class T>
inline void ColMajorInd2Sub(const size_t& n_rows, const T<size_t>& ind, T<size_t>& row, T<size_t>& col)
{
	RowMajorInd2Sub(n_rows, ind, col, row);
}






template <class T>
Matrix<T> ColMajor2RowMajor(const size_t& m, const size_t& n, const T* const &p_ColMajor)
{
	Matrix<T> Res(m, n);

	for (size_t i=0; i<m;i++)
		for (size_t j=0; j<n;j++)
			Res[i][j]=p_ColMajor[j*m+i];

	return Res;
}

template <class T>
Matrix<complex<T> > ComplexColMajor2RowMajor(const size_t& m, const size_t& n, const T* const &p_ColMajor_real, const T* const &p_ColMajor_imaginary)
{
	Matrix<complex<T> > Res(m, n);

	for (size_t i=0; i<m;i++)
		for (size_t j=0; j<n;j++)
			Res[i][j]=complex<T>(p_ColMajor_real[j*m+i],p_ColMajor_imaginary[j*m+i]);

	return Res;
}


template <class T>
void RowMajor2ColMajor(const Matrix<T>& H_RowMajor, T* const &p_ColMajor)
{
	for (size_t j=0; j<H_RowMajor.cols;j++)
		for (size_t i=0; i<H_RowMajor.rows;i++)
			p_ColMajor[j*H_RowMajor.rows+i]=H_RowMajor[i][j];
}

template <class T>
void ComplexRowMajor2ColMajor(const Matrix<complex<T> >& H_RowMajor, T* const &p_ColMajor_real,  T* const &p_ColMajor_imaginary)
{
	for (size_t j=0; j<H_RowMajor.cols;j++)
		for (size_t i=0; i<H_RowMajor.rows;i++)
		{
			p_ColMajor_real[j*H_RowMajor.rows+i]=H_RowMajor[i][j].real();
			p_ColMajor_imaginary[j*H_RowMajor.rows+i]=H_RowMajor[i][j].imag();
		}
}


#endif