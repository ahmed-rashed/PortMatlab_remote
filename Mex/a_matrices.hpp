#ifndef a_matrices_h
#define a_matrices_h

#include "a_matrices_internals.hpp"
#include <valarray>
#include <complex>
#include "valarray_ComplexRealOperators.hpp"
#include <iostream>

namespace a_matrices{

template <class T>
class Matrix
{
private:
	std::valarray<T> va;

protected:
	//This constructor is not allowed
	Matrix() {}

public:
	const size_t rows, cols;

	//Creates a Matrix with rrows*ccols elements valarray default constructor
	Matrix(const size_t rrows, const size_t ccols)
		:rows(rrows), cols(ccols), va(std::valarray<T> (rrows*ccols)) {}

	//Creates a Matrix with rrows*ccols elements valarray initialized by value.
	Matrix (const size_t rrows, const size_t ccols, const T& value)
		:rows(rrows), cols(ccols), va(std::valarray<T> (value, rrows*ccols)) {}

	//Creates a Matrix with rrows*ccols elements valarray initialized by the values of the elements in *parray.
	//The caller must ensure that *parray contains rrows*ccols elements; otherwise, the behavior is undefined.
	Matrix (const size_t rrows, const size_t ccols, const T* parray)
		:rows(rrows), cols(ccols), va(std::valarray<T> (parray, rrows*ccols)) {}

	//Creates a Matrix with rrows*ccols elements valarray initialized by another valarray.
	Matrix (const size_t rrows, const size_t ccols, const std::valarray<T> valarr)
		:rows(rrows), cols(ccols)
	{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
		if (valarr.size()!=(rows*cols)) throw("Valarr must have size=rows*cols.");
#endif
		va=valarr;
	}

	//The copy constructor is not needed

	//destructor not needed

	std::valarray<T> va_copy() {return va;}

	size_t size() const	{return va.size();}

	//assignment operator
	Matrix<T>& operator=(const Matrix<T>& Mat)
	{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
		if (!equal_size(*this,Mat)) throw("For matrix assignment, matrices dimension must agree.");
#endif
		va=Mat.va;
		return *this;
	}

	Matrix<T>& operator=(const T& value)
	{
		va=value;
		return *this;
	}

	Matrix<T>& operator=(const std::gslice_array<T>& Mat)
	{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
		if (!equal_size(*this,Mat)) throw("For matrix assignment, matrices dimension must agree.");
#endif
		va=Mat.va;
		return *this;
	}

	//Returns pointer to row i, i.e., 1D array and not Matrix
	//When executing Matrix[i][j] the following function overloads only the 1st operator [i]
	//The second [] operator is not overloaded as it affects the result of this function
	//which is a pointer to row i, i.e., 1D array and not Matrix
	//The returned pointer not declared precisely. Its precise declaration is T(*)[cols]
	//This assumes the exsistance of second subscribt, [j], in the calling command.
	//So this function, Matrix[i], returns pointer and not reference.
	//The default [] operator, (Matrix[i])[j], then return the correct result suitable for l-value or r-value
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
	const next_sub<T> operator[](const size_t &i) const
#else
	const T* operator[](const size_t &i) const
#endif
	{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
		if (i>=rows) throw("Row subscript is out of bounds.");
		return next_sub<T>(&(((const_cast<Matrix<T>*>(this))->va)[i*cols]),cols);
//		return next_sub<T>(&(va[i*cols]),cols);//Didn't work
#else
		return &(((const_cast<Matrix<T>*>(this))->va)[i*cols]);
//		return &(va[i*cols]);//Didn't work
#endif
	}

#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
	next_sub<T> operator[](const size_t &i)
#else
	T* operator[](const size_t &i)
#endif
	{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
		if (i>=rows) throw("Row subscript is out of bounds.");
		return next_sub<T>(&(va[i*cols]),cols);
#else
		return &(va[i*cols]);
#endif
	}

	//Vector of indices
	std::indirect_array<T> operator[](const std::valarray<size_t> &Ind)
	{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
		if (any(Ind>=size())) throw("Subscripts must be within matrix bounds.");
#endif
		return va[Ind];
	}

	//Vector of indices
	std::valarray<T> operator[](const std::valarray<size_t> &Ind) const
	{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
		if (any(Ind>=size())) throw("Subscripts must be within matrix bounds.");
#endif
		return va[Ind];
//		return ((const_cast<Matrix<T>*>(this))->va)[Ind];
	}

	//Logigal indexing
	std::mask_array<T> operator[](const Matrix<bool> &Msk)
	{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
		if (!equal_size(*this,Msk)) throw("Mask array must has the same dimensions as the masked array.");
#endif
		return va[Msk.va];
	}

	//Logigal indexing
	std::valarray<T> operator[](const std::valarray<T> &Msk) const
	{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
		if (!equal_size(*this,Msk)) throw("Mask array must has the same dimensions as the masked array.");
#endif
		return va[Msk.va];
//		return ((const_cast<Matrix<T>*>(this))->va)[Msk.va];
	}


	//row
	std::slice_array<T> row(size_t i)
	{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
		if (i>=rows) throw("Row subscript is out of bounds.");
#endif		
		return va[std::slice (i*cols,cols,1)];
	}

	//row
	std::valarray<T> row(size_t i) const
	{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
		if (i>=rows) throw("Row subscript is out of bounds.");
#endif		
		return va[std::slice (i*cols,cols,1)];
	}

	//col
	std::slice_array<T> col(size_t j)
	{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
		if (j>=cols) throw("Column subscript is out of bounds.");
#endif		
		return va[std::slice (j,rows,cols)];
	}

	//col
	std::valarray<T> col(size_t j) const
	{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
		if (j>=cols) throw("Column subscript is out of bounds.");
#endif		
		return va[std::slice (j,rows,cols)];
	}

	//SubMatrix
	std::gslice_array<T> SubMatrix(size_t i_start, size_t i_end, size_t j_start, size_t j_end)
	{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
		if ((i_start>=rows) || (i_end>=rows) || (j_start>=cols) || (j_end>=cols)) throw("Indices must be within matrix bounds.");
#endif
		size_t len[]={i_end-i_start+1, j_end-j_start+1};
		size_t str[]={cols, 1};
		std::valarray<size_t> lengths(len,2);
		std::valarray<size_t> strids(str,2);
		std::gslice gs(i_start*cols+j_start,lengths,strids);
		
		return va[gs];
	}

	//SubMatrix
	std::valarray<T> SubMatrix(size_t i_start, size_t i_end, size_t j_start, size_t j_end) const
	{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
		if ((i_start>=rows) || (i_end>=rows) || (j_start>=cols) || (j_end>=cols)) throw("Indices must be within matrix bounds.");
#endif
		size_t len[]={i_end-i_start+1, j_end-j_start+1};
		size_t str[]={cols, 1};
		std::valarray<size_t> lengths(len,2);
		std::valarray<size_t> strids(str,2);
		std::gslice gs(i_start*cols+j_start,lengths,strids);
		
		return va[gs];
	}

	//SubMatrix
	std::gslice_array<T> SubMatrix(size_t i_start, size_t i_end, size_t i_strid, size_t j_start, size_t j_end, size_t j_strid)
	{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
		if ((i_start>=rows) || (i_end>=rows) || (j_start>=cols) || (j_end>=cols)) throw("Indices must be within matrix bounds.");
#endif
		size_t len[]={floor(static_cast<double>(i_end-i_start)/i_strid)+1, floor(static_cast<double>(j_end-j_start)/j_strid)+1};
		size_t str[]={i_strid*cols, j_strid};
		std::valarray<size_t> lengths(len,2);
		std::valarray<size_t> strids(str,2);
		std::gslice gs(i_start*cols+j_start,lengths,strids);
		
		return va[gs];
	}

	//SubMatrix
	std::valarray<T> SubMatrix(size_t i_start, size_t i_end, size_t i_strid, size_t j_start, size_t j_end, size_t j_strid) const
	{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
		if ((i_start>=rows) || (i_end>=rows) || (j_start>=cols) || (j_end>=cols)) throw("Indices must be within matrix bounds.");
#endif
		{floor(static_cast<double>(i_end-i_start)/i_strid)+1, floor(static_cast<double>(j_end-j_start)/j_strid)+1}
		size_t str[]={i_strid*cols, j_strid};
		std::valarray<size_t> lengths(len,2);
		std::valarray<size_t> strids(str,2);
		std::gslice gs(i_start*cols+j_start,lengths,strids);
		
		return va[gs];
	}

	//unary operators
	Matrix<T> operator-() const
	{
		return Matrix<T>(rows,cols,-va);
	}

	Matrix<T> operator!() const
	{
		return Matrix<T>(rows,cols,!va);
	}

 	template <class T> friend Matrix<T> operator+(const Matrix<T> &, const Matrix<T> &);
	template <class T> friend Matrix<T> operator+(const Matrix<T> &, const T &);
 	template <class T> friend Matrix<T> operator-(const Matrix<T> &, const Matrix<T> &);
	template <class T> friend Matrix<T> operator-(const Matrix<T> &, const T &);
	template <class T> friend Matrix<T> operator-(const T &, const Matrix<T> &);
 	template <class T> friend Matrix<T> operator*(const Matrix<T> &, const Matrix<T> &);
	template <class T> friend Matrix<T> operator*(const Matrix<T> &, const T &);
 	template <class T> friend Matrix<T> operator/(const Matrix<T> &, const Matrix<T> &);
	template <class T> friend Matrix<T> operator/(const Matrix<T> &, const T &);
	template <class T> friend Matrix<T> operator/(const T &, const Matrix<T> &);


 	template <class T> friend Matrix<std::complex<T> > operator+(const Matrix<std::complex<T> > &, const Matrix<T> &);
	template <class T> friend Matrix<std::complex<T> > operator+(const Matrix<T> &, const Matrix<std::complex<T> > &);
	template <class T> friend Matrix<std::complex<T> > operator+(const Matrix<std::complex<T> > &, const T &);
	template <class T> friend Matrix<std::complex<T> > operator+(const Matrix<T> &, const std::complex<T> &);
 	template <class T> friend Matrix<std::complex<T> > operator-(const Matrix<std::complex<T> > &, const Matrix<T> &);
	template <class T> friend Matrix<std::complex<T> > operator-(const Matrix<T> &, const Matrix<std::complex<T> > &);
	template <class T> friend Matrix<std::complex<T> > operator-(const Matrix<std::complex<T> > &, const T &);
	template <class T> friend Matrix<std::complex<T> > operator-(const Matrix<T> &, const std::complex<T> &);
	template <class T> friend Matrix<std::complex<T> > operator-(const std::complex<T> &, const Matrix<T> &);
	template <class T> friend Matrix<std::complex<T> > operator-(const T &, const Matrix<std::complex<T> > &);
 	template <class T> friend Matrix<std::complex<T> > operator*(const Matrix<std::complex<T> > &, const Matrix<T> &);
	template <class T> friend Matrix<std::complex<T> > operator*(const Matrix<T> &, const Matrix<std::complex<T> > &);
	template <class T> friend Matrix<std::complex<T> > operator*(const Matrix<std::complex<T> > &, const T &);
	template <class T> friend Matrix<std::complex<T> > operator*(const Matrix<T> &, const std::complex<T> &);
 	template <class T> friend Matrix<std::complex<T> > operator/(const Matrix<std::complex<T> > &, const Matrix<T> &);
	template <class T> friend Matrix<std::complex<T> > operator/(const Matrix<T> &, const Matrix<std::complex<T> > &);
	template <class T> friend Matrix<std::complex<T> > operator/(const Matrix<std::complex<T> > &, const T &);
	template <class T> friend Matrix<std::complex<T> > operator/(const Matrix<T> &, const std::complex<T> &);
	template <class T> friend Matrix<std::complex<T> > operator/(const std::complex<T> &, const Matrix<T> &);
	template <class T> friend Matrix<std::complex<T> > operator/(const T &, const Matrix<std::complex<T> > &);
	template <class T> friend Matrix<T> real(const Matrix<std::complex<T> >& complexMat);
	template <class T> friend Matrix<T> imag(const Matrix<std::complex<T> >& complexMat);
	template <class T> friend Matrix<std::complex<T> > conj(const Matrix<std::complex<T> >& complexMat);
	template <class T> friend Matrix<T> abs(const Matrix<std::complex<T> >& complexMat);
	template <class T> friend Matrix<T> arg(const Matrix<std::complex<T> >& complexMat);


	template <class T> friend Matrix<bool> operator==(const Matrix<T>&, const Matrix<T>&);
	template <class T> friend Matrix<bool> operator==(const Matrix<T>&, const T&);
	template <class T> friend Matrix<bool> operator!=(const Matrix<T>&, const Matrix<T>&);
	template <class T> friend Matrix<bool> operator!=(const Matrix<T>&, const T&);
	template <class T> friend Matrix<bool> operator<(const Matrix<T>&, const Matrix<T>&);
	template <class T> friend Matrix<bool> operator<(const Matrix<T>&, const T&);
	template <class T> friend Matrix<bool> operator<=(const Matrix<T>&, const Matrix<T>&);
	template <class T> friend Matrix<bool> operator<=(const Matrix<T>&, const T&);
	template <class T> friend Matrix<bool> operator>(const Matrix<T>&, const Matrix<T>&);
	template <class T> friend Matrix<bool> operator>(const Matrix<T>&, const T&);
	template <class T> friend Matrix<bool> operator>=(const Matrix<T>&, const Matrix<T>&);
	template <class T> friend Matrix<bool> operator>=(const Matrix<T>&, const T&);
	template <class T> friend Matrix<bool> operator&&(const Matrix<T>&, const Matrix<T>&);
	template <class T> friend Matrix<bool> operator&&(const Matrix<T>&, const T&);
	template <class T> friend Matrix<bool> operator||(const Matrix<T>&, const Matrix<T>&);
	template <class T> friend Matrix<bool> operator||(const Matrix<T>&, const T&);
	template <class T> friend Matrix<T> abs(const Matrix<T>& );
	template <class T> friend Matrix<T> sqrt(const Matrix<T>& );
	template <class T> friend Matrix<T> pow(const Matrix<T>&, const Matrix<T>&);
	template <class T> friend Matrix<T> pow(const Matrix<T>&, const T&);
	//template <class T> friend T sum(const Matrix<T>& );
	template <class T> friend Matrix<T> scale_rows(const Matrix<T>& ,const Matrix<T>& );
	template <class T> friend Matrix<T> scale_cols(const Matrix<T>& ,const Matrix<T>& );
	friend bool any(const Matrix<bool> & );
	friend bool all(const Matrix<bool> & );
	friend inline std::valarray<size_t> find(Matrix<bool> );
	template <class T> friend std::valarray<T> MatDiag(const Matrix<T>&);
	template <class T> friend Matrix<T> DiagMat(const std::valarray<T> &);

}; //End of the Matrix class


///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

//Addition operator
template <class T>
inline Matrix<T> operator+(const Matrix<T> &Mat1, const Matrix<T> &Mat2)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
		if (!equal_size(Mat1,Mat2)) throw("Matrices must have identical size for addition.");
#endif
	return Matrix<T>(Mat1.rows,Mat1.cols,Mat1.va+Mat2.va);
}

template <class T>
inline Matrix<T> operator+(const Matrix<T> &Mat, const T &value)
{
	return Matrix<T>(Mat.rows,Mat.cols,Mat.va+value);
}

template <class T>
inline Matrix<T> operator+(const T &value, const Matrix<T> &Mat) {return Mat+value;}

//Subtraction operator
template <class T>
inline Matrix<T> operator-(const Matrix<T> &Mat1, const Matrix<T> &Mat2)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
		if (!equal_size(Mat1,Mat2)) throw("Matrices must have identical size for subtraction.");
#endif
	return Matrix<T>(Mat1.rows,Mat1.cols,Mat1.va-Mat2.va);
}

template <class T>
inline Matrix<T> operator-(const Matrix<T> &Mat, const T &value)
{
	return Matrix<T>(Mat.rows,Mat.cols,Mat.va-value);
}

template <class T>
inline Matrix<T> operator-(const T &value, const Matrix<T> &Mat)
{
	return Matrix<T>(Mat.rows,Mat.cols,value-Mat.va);
}

//Multiplication operator
template <class T>
inline Matrix<T> operator*(const Matrix<T> &Mat1, const Matrix<T> &Mat2)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
	if (!equal_size(Mat1,Mat2)) throw("Matrices must have identical size for multiplication.");
#endif
	return Matrix<T>(Mat1.rows,Mat1.cols,Mat1.va*Mat2.va);
}

template <class T>
inline Matrix<T> operator*(const Matrix<T> &Mat, const T &value)
{
	return Matrix<T>(Mat.rows,Mat.cols,Mat.va*value);
}

template <class T>
inline Matrix<T> operator*(const T &value, const Matrix<T> &Mat) {return Mat*value;}

//Division operator
template <class T>
inline Matrix<T> operator/(const Matrix<T> &Mat1, const Matrix<T> &Mat2)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
	if (!equal_size(Mat1,Mat2)) throw("Matrices must have identical size for division.");
#endif
	return Matrix<T>(Mat1.rows,Mat1.cols,Mat1.va/Mat2.va);
}

template <class T>
inline Matrix<T> operator/(const Matrix<T> &Mat, const T &value)
{
	return Matrix<T>(Mat.rows,Mat.cols,Mat.va/value);
}

template <class T>
inline Matrix<T> operator/(const T &value, const Matrix<T> &Mat)
{
	return Matrix<T>(Mat.rows,Mat.cols,value/Mat.va);
}

//logical operators
//==
template <class T>
inline Matrix<bool> operator== (const Matrix<T>& Mat1, const Matrix<T>& Mat2)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
	if (!equal_size(Mat1,Mat2)) throw("Matrices must have identical size for ==.");
#endif
	return Matrix<bool>(Mat1.rows,Mat1.cols,Mat1.va==Mat2.va);
}

template <class T>
inline Matrix<bool> operator== (const Matrix<T>& Mat, const T& value)
{
	return Matrix<bool>(Mat.rows,Mat.cols,Mat.va==value);
}

template <class T>
inline Matrix<bool> operator==(const T &value, const Matrix<T> &Mat) {return Mat==value;}

//!=
template <class T>
inline Matrix<bool> operator!= (const Matrix<T>& Mat1, const Matrix<T>& Mat2)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
	if (!equal_size(Mat1,Mat2)) throw("Matrices must have identical size for !=.");
#endif
	return Matrix<bool>(Mat1.rows,Mat1.cols,Mat1.va!=Mat2.va);
}

template <class T>
inline Matrix<bool> operator!= (const Matrix<T>& Mat, const T& value)
{
	return Matrix<bool>(Mat.rows,Mat.cols,Mat.va!=value);
}

template <class T>
inline Matrix<bool> operator!=(const T &value, const Matrix<T> &Mat) {return Mat!=value;}

//<
template <class T>
inline Matrix<bool> operator< (const Matrix<T>& Mat1, const Matrix<T>& Mat2)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
	if (!equal_size(Mat1,Mat2)) throw("Matrices must have identical size for <.");
#endif
	return Matrix<bool>(Mat1.rows,Mat1.cols,Mat1.va<Mat2.va);
}

template <class T>
inline Matrix<bool> operator< (const Matrix<T>& Mat, const T& value)
{
	return Matrix<bool>(Mat.rows,Mat.cols,Mat.va<value);
}

template <class T>
inline Matrix<bool> operator<(const T &value, const Matrix<T> &Mat) {return Mat<value;}

//<=
template <class T>
inline Matrix<bool> operator<= (const Matrix<T>& Mat1, const Matrix<T>& Mat2)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
	if (!equal_size(Mat1,Mat2)) throw("Matrices must have identical size for <=.");
#endif
	return Matrix<bool>(Mat1.rows,Mat1.cols,Mat1.va<=Mat2.va);
}

template <class T>
inline Matrix<bool> operator<= (const Matrix<T>& Mat, const T& value)
{
	return Matrix<bool>(Mat.rows,Mat.cols,Mat.va<=value);
}

template <class T>
inline Matrix<bool> operator<=(const T &value, const Matrix<T> &Mat) {return Mat<=value;}

//>
template <class T>
inline Matrix<bool> operator> (const Matrix<T>& Mat1, const Matrix<T>& Mat2)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
	if (!equal_size(Mat1,Mat2)) throw("Matrices must have identical size for >.");
#endif
	return Matrix<bool>(Mat1.rows,Mat1.cols,Mat1.va>Mat2.va);
}

template <class T>
inline Matrix<bool> operator> (const Matrix<T>& Mat, const T& value)
{
	return Matrix<bool>(Mat.rows,Mat.cols,Mat.va>value);
}

template <class T>
inline Matrix<bool> operator>(const T &value, const Matrix<T> &Mat) {return Mat>value;}

//>=
template <class T>
inline Matrix<bool> operator>= (const Matrix<T>& Mat1, const Matrix<T>& Mat2)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
	if (!equal_size(Mat1,Mat2)) throw("Matrices must have identical size for >=.");
#endif
	return Matrix<bool>(Mat1.rows,Mat1.cols,Mat1.va>=Mat2.va);
}

template <class T>
inline Matrix<bool> operator>= (const Matrix<T>& Mat, const T& value)
{
	return Matrix<bool>(Mat.rows,Mat.cols,Mat.va>=value);
}

template <class T>
inline Matrix<bool> operator>=(const T &value, const Matrix<T> &Mat) {return Mat>=value;}

//||
template <class T>
inline Matrix<bool> operator|| (const Matrix<T>& Mat1, const Matrix<T>& Mat2)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
	if (!equal_size(Mat1,Mat2)) throw("Matrices must have identical size for ||.");
#endif
	return Matrix<bool>(Mat1.rows,Mat1.cols,Mat1.va||Mat2.va);
}

template <class T>
inline Matrix<bool> operator|| (const Matrix<T>& Mat, const T& value)
{
	return Matrix<bool>(Mat.rows,Mat.cols,Mat.va||value);
}

template <class T>
inline Matrix<bool> operator||(const T &value, const Matrix<T> &Mat) {return Mat||value;}

//&&
template <class T>
inline Matrix<bool> operator&& (const Matrix<T>& Mat1, const Matrix<T>& Mat2)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
	if (!equal_size(Mat1,Mat2)) throw("Matrices must have identical size for &&.");
#endif
	return Matrix<bool>(Mat1.rows,Mat1.cols,Mat1.va&&Mat2.va);
}

template <class T>
inline Matrix<bool> operator&& (const Matrix<T>& Mat, const T& value)
{
	return Matrix<bool>(Mat.rows,Mat.cols,Mat.va&&value);
}

template <class T>
inline Matrix<bool> operator&&(const T &value, const Matrix<T> &Mat) {return Mat&&value;}

template <class T>
Matrix<T> abs(const Matrix<T>& Mat)
{
	return Matrix<T>(Mat.rows, Mat.cols,abs(Mat.va));
}

template <class T>
Matrix<T> sqrt(const Matrix<T>& Mat)
{
	return Matrix<T>(Mat.rows, Mat.cols,sqrt(Mat.va));
}

template <class T>
Matrix<T> pow(const Matrix<T>& Mat1, const Matrix<T>&  Mat2)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
	if (!equal_size(Mat1,Mat2)) throw("Matrices must have identical size for &&.");
#endif
	return Matrix<T>(Mat1.rows, Mat1.cols,pow(Mat1.va, Mat2.va));
}

template <class T>
Matrix<T> pow(const Matrix<T>& Mat, const T& val)
{
	return Matrix<T>(Mat.rows, Mat.cols,pow(Mat.va, val));
}

template <class T>
Matrix<T> pow(const T& val, const Matrix<T>& Mat)
{
	return pow(Mat, val);
}

/*template <class T>
T sum(const Matrix<T>& Mat)
{
	if ((Mat.rows!=1) && (Mat.cols!=1)) throw("Input must be row vector or column vector.");
	return Mat.va.sum();
}*/

template <class T>
Matrix<T> scale_rows(const Matrix<T>& Mat,const std::valarray<T>& Scale_vec)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
	if (Mat.rows!=Scale_vec.size()) throw("Scale_vec must have same number of elements as rows of Mat.");
#endif
	Matrix<T> Res(Mat.rows, Mat.cols);
	for (size_t i=0;i<Mat.rows;i++)
		Res.SubMatrix(i,i,0,Mat.cols-1)=static_cast<std::valarray<T> >(Mat.SubMatrix(i,i,0,Mat.cols-1))*Scale_vec[i];
	
	return Res;
}

template <class T>
Matrix<T> scale_rows(const Matrix<T>& Mat,const Matrix<T>& Scale_col)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
	if ((Scale_col.rows != Mat.rows) || (Scale_col.cols != 1)) throw("Scale_col must be column vector with same number of rows as Mat.");
#endif	
	return scale_rows( Mat, Scale_col.va);
}

template <class T>
Matrix<T> scale_cols(const Matrix<T>& Mat,const std::valarray<T>& Scale_vec)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
	if (Mat.cols!=Scale_vec.size()) throw("Scale_vec must have same number of elements as columns of Mat.");
#endif
	Matrix<T> Res(Mat.rows, Mat.cols);
	for (size_t i=0;i<Mat.cols;i++)
		Res.SubMatrix(0,Mat.rows-1,i,i)=static_cast<std::valarray<T> >(Mat.SubMatrix(0,Mat.cols-1,i,i))*Scale_vec[i];
	
	return Res;
}

template <class T>
Matrix<T> scale_cols(const Matrix<T>& Mat,const Matrix<T>& Scale_col)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
	if ((Scale_col.rows != Mat.cols) || (Scale_col.cols != 1)) throw("Scale_col must be column vector with same number of rows as Mat columns.");
#endif	
	return scale_cols(Mat, Scale_col.va);
}

template <class T>
Matrix<T> transpose(const Matrix<T> Mat)
{
	Matrix<T> Res(Mat.cols, Mat.rows);
	for (size_t i=0;i<Res.rows;i++)
		for (size_t j=0;j<Res.cols;j++)
			Res[i][j]=Mat[j][i];

	return Res;

	/*	
	size_t len[]={Mat.cols, Mat.rows};
	size_t str[]={1, Mat.cols};
	std::valarray<size_t> lengths(len,2);
	std::valarray<size_t> strids(str,2);
	std::gslice gs(0,lengths,strids);
	Res.va=static_cast<std::valarray<T> >(Mat.va[gs]);
	*/
}

inline bool any(const std::valarray<bool> &Carr)
{
	for(size_t i = 0;i<Carr.size();i++)
		if (Carr[i]) return true;
	return false;
}

inline bool any(const Matrix<bool> &Mat) {return any(Mat.va);}

inline bool all(const std::valarray<bool> &Carr)
{
	for(size_t i = 0;i<Carr.size();i++)
		if (!Carr[i]) return false;
	return true;
}

inline bool all(const Matrix<bool> &Mat) {return all(Mat.va);}

template <class T1, class T2>
inline bool equal_size(const Matrix<T1> &a, const Matrix<T2> &b)
{
	if ((a.rows==b.rows) && (a.cols==b.cols))
		return true;
	else
		return false;
}

inline std::valarray<size_t> find(std::valarray<bool> Arr)
{
	size_t counter=0;
	for (size_t i=0;i<Arr.size();i++) if (Arr[i]) counter++;

	std::valarray<size_t> indx(counter);
	
	counter=0;
	for (size_t i=0;i<Arr.size();i++){
		if (Arr[i]){
			indx[counter]=i;
			counter++;
		}
	}
#if defined(_DEBUG)
		if (counter!=indx.size()) throw("This function don't run correctly.");
#endif
		return indx;
}

inline std::valarray<size_t> find(Matrix<bool> Mat){return find(Mat.va);}

/*template <class T>
Matrix<T> sort(const Matrix<T>& Mat, const enum {1, 2} dim, Matrix<size_t> Indx)
{
	Matrix<T> Res(Mat.rows, Mat.cols);
	if (dim==1) //sort each column
	{
		std::valarray<T> Val_tmp(Mat.rows);
		for (size_t i=0; i<Mat.cols; i++)
		{
			val_tmp=Mat.col(i)
		}
	}

}*/


/*template <class T>
void printMatrix (const Matrix<T>& va)
{
	std::cout << std::endl << "Matrix:" << std::endl;
    for (size_t i=0; i<va.rows; i++) {
		for (size_t j=0; j<va.cols; j++){
			T var=va[i][j];
			std::cout << var << ' ';
		}
		std::cout << std::endl;
    }
    std::cout << std::endl;
}
*/

template <class T>
std::ostream & operator << (std::ostream& stream, Matrix<T>& Mat)
{

	stream << std::endl;
    for (size_t i=0; i<Mat.rows; i++)
	{
		stream << Mat[i][0];
		for (size_t j=1; j<Mat.cols; j++) stream << ' ' << Mat[i][j];

		if (i < Mat.rows-1) stream << std::endl;
    }

	return stream;
}





///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Addition operator
template <class T>
inline Matrix<std::complex<T> > operator+(const Matrix<std::complex<T> > &Mat1, const Matrix<T> &Mat2)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
		if (!equal_size(Mat1,Mat2)) throw("Matrices must have identical size for addition.");
#endif
	return Matrix<std::complex<T> >(Mat1.rows,Mat1.cols,Mat1.va+Mat2.va);
}

template <class T>
inline Matrix<std::complex<T> > operator+(const Matrix<T> &Mat1, const Matrix<std::complex<T> > &Mat2)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
		if (!equal_size(Mat1,Mat2)) throw("Matrices must have identical size for addition.");
#endif
	return Matrix<std::complex<T> >(Mat1.rows,Mat1.cols,Mat1.va+Mat2.va);
}

template <class T>
inline Matrix<std::complex<T> > operator+(const Matrix<std::complex<T> > &Mat, const T &value)
{
	return Matrix<std::complex<T> >(Mat.rows,Mat.cols,Mat.va+value);
}

template <class T>
inline Matrix<std::complex<T> > operator+(const Matrix<T> &Mat, const std::complex<T> &value)
{
	return Matrix<std::complex<T> >(Mat.rows,Mat.cols,Mat.va+value);
}

template <class T>
inline Matrix<std::complex<T> > operator+(const std::complex<T> &value, const Matrix<T> &Mat) {return Mat+value;}

template <class T>
inline Matrix<std::complex<T> > operator+(const T &value, const Matrix<std::complex<T> > &Mat) {return Mat+value;}

//Subtraction operator
template <class T>
inline Matrix<std::complex<T> > operator-(const Matrix<std::complex<T> > &Mat1, const Matrix<T> &Mat2)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
		if (!equal_size(Mat1,Mat2)) throw("Matrices must have identical size for subtraction.");
#endif
	return Matrix<std::complex<T> >(Mat1.rows,Mat1.cols,Mat1.va-Mat2.va);
}

template <class T>
inline Matrix<std::complex<T> > operator-(const Matrix<T> &Mat1, const Matrix<std::complex<T> > &Mat2)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
		if (!equal_size(Mat1,Mat2)) throw("Matrices must have identical size for subtraction.");
#endif
	return Matrix<std::complex<T> >(Mat1.rows,Mat1.cols,Mat1.va-Mat2.va);
}

template <class T>
inline Matrix<std::complex<T> > operator-(const Matrix<std::complex<T> > &Mat, const T &value)
{
	return Matrix<std::complex<T> >(Mat.rows,Mat.cols,Mat.va-value);
}

template <class T>
inline Matrix<std::complex<T> > operator-(const Matrix<T> &Mat, const std::complex<T> &value)
{
	return Matrix<std::complex<T> >(Mat.rows,Mat.cols,Mat.va-value);
}

template <class T>
inline Matrix<std::complex<T> > operator-(const std::complex<T> &value, const Matrix<T> &Mat)
{
	return Matrix<std::complex<T> >(Mat.rows,Mat.cols,value-Mat.va);
}

template <class T>
inline Matrix<std::complex<T> > operator-(const T &value, const Matrix<std::complex<T> > &Mat)
{
	return Matrix<std::complex<T> >(Mat.rows,Mat.cols,value-Mat.va);
}

//Multiplication operator
template <class T>
inline Matrix<std::complex<T> > operator*(const Matrix<std::complex<T> > &Mat1, const Matrix<T> &Mat2)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
	if (!equal_size(Mat1,Mat2)) throw("Matrices must have identical size for multiplication.");
#endif
	return Matrix<std::complex<T> >(Mat1.rows,Mat1.cols,Mat1.va*Mat2.va);
}

template <class T>
inline Matrix<std::complex<T> > operator*(const Matrix<T> &Mat1, const Matrix<std::complex<T> > &Mat2)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
	if (!equal_size(Mat1,Mat2)) throw("Matrices must have identical size for multiplication.");
#endif
	return Matrix<std::complex<T> >(Mat1.rows,Mat1.cols,Mat1.va*Mat2.va);
}

template <class T>
inline Matrix<std::complex<T> > operator*(const Matrix<std::complex<T> > &Mat, const T &value)
{
	return Matrix<std::complex<T> >(Mat.rows,Mat.cols,Mat.va*value);
}

template <class T>
inline Matrix<std::complex<T> > operator*(const Matrix<T> &Mat, const std::complex<T> &value)
{
	return Matrix<std::complex<T> >(Mat.rows,Mat.cols,Mat.va*value);
}

template <class T>
inline Matrix<std::complex<T> > operator*(const std::complex<T> &value, const Matrix<T> &Mat) {return Mat*value;}

template <class T>
inline Matrix<std::complex<T> > operator*(const T &value, const Matrix<std::complex<T> > &Mat) {return Mat*value;}

//Division operator
template <class T>
inline Matrix<std::complex<T> > operator/(const Matrix<std::complex<T> > &Mat1, const Matrix<T> &Mat2)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
	if (!equal_size(Mat1,Mat2)) throw("Matrices must have identical size for division.");
#endif
	return Matrix<std::complex<T> >(Mat1.rows,Mat1.cols,Mat1.va/Mat2.va);
}

template <class T>
inline Matrix<std::complex<T> > operator/(const Matrix<T> &Mat1, const Matrix<std::complex<T> > &Mat2)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
	if (!equal_size(Mat1,Mat2)) throw("Matrices must have identical size for division.");
#endif
	return Matrix<std::complex<T> >(Mat1.rows,Mat1.cols,Mat1.va/Mat2.va);
}

template <class T>
inline Matrix<std::complex<T> > operator/(const Matrix<std::complex<T> > &Mat, const T &value)
{
	return Matrix<std::complex<T> >(Mat.rows,Mat.cols,Mat.va/value);
}

template <class T>
inline Matrix<std::complex<T> > operator/(const Matrix<T> &Mat, const std::complex<T> &value)
{
	return Matrix<std::complex<T> >(Mat.rows,Mat.cols,Mat.va/value);
}

template <class T>
inline Matrix<std::complex<T> > operator/(const std::complex<T> &value, const Matrix<T> &Mat)
{
	return Matrix<std::complex<T> >(Mat.rows,Mat.cols,value/Mat.va);
}

template <class T>
inline Matrix<std::complex<T> > operator/(const T &value, const Matrix<std::complex<T> > &Mat)
{
	return Matrix<std::complex<T> >(Mat.rows,Mat.cols,value/Mat.va);
}

template <class T>
Matrix<T> real(const Matrix<std::complex<T> >& complexMat)
{
	return Matrix<T>(complexMat.rows,complexMat.cols,real(complexMat.va));
}

template <class T>
Matrix<T> imag(const Matrix<std::complex<T> >& complexMat)
{
	return Matrix<T>(complexMat.rows,complexMat.cols,imag(complexMat.va));
}

template <class T>
Matrix<std::complex<T> > conj(const Matrix<std::complex<T> >& complexMat)
{
	return Matrix<std::complex<T> >(complexMat.rows,complexMat.cols,conj(complexMat.va));
}

template <class T>
Matrix<T> abs(const Matrix<std::complex<T> >& complexMat)
{
	return Matrix<T>(complexMat.rows,complexMat.cols,abs(complexMat.va));
}

template <class T>
Matrix<T> arg(const Matrix<std::complex<T> >& complexMat)
{
	return Matrix<T>(complexMat.rows,complexMat.cols,arg(complexMat.va));
}

template <class T>
std::valarray<T> MatDiag(const Matrix<T> &Mat)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
		if (Mat.rows!=Mat.cols) throw("Matrix must be square.");
#endif
	return Mat.va[std::slice(0,Mat.rows,Mat.rows+1)];
}

template <class T>
Matrix<T> DiagMat(const std::valarray<T> &vec)
{
	size_t N=vec.size();
	Matrix<T> Mat(N,N);
	Mat.va[std::slice(0,N,N+1)]=vec;

	return Mat;
}

}//end of the namespace a_matrices
#endif