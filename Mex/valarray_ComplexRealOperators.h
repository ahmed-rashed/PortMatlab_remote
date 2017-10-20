#ifndef valarray_ComplexRealOperatore_h
#define valarray_ComplexRealOperatore_h

#include <complex>
#include <valarray>

//Addition operator
template <class T>
inline std::valarray<std::complex<T> > operator+(const std::valarray<std::complex<T> > &varray1, const std::valarray<T> &varray2)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
		if (varray1.size()!=varray2.size()) throw("operator+ valarray operands must have the same size.");
#endif

	std::valarray<std::complex<T> > Result(varray1);
	for (int n=0; n<varray1.size();n++)
		Result[n]+=varray2[n];

	return Result;
}

template <class T>
inline std::valarray<std::complex<T> > operator+(const std::valarray<T> &varray2, const std::valarray<std::complex<T> > &varray1) {return varray1+varray2;}

template <class T>
inline std::valarray<std::complex<T> > operator+(const std::valarray<std::complex<T> > &varray, const T &value)
{
	std::valarray<std::complex<T> > Result(varray);
	Result+=std::complex<T>(value,0);

	return Result;
}

template <class T>
inline std::valarray<std::complex<T> > operator+(const T &value, const std::valarray<std::complex<T> > &varray) {return varray+value;}

template <class T>
inline std::valarray<std::complex<T> > operator+(const std::valarray<T> &varray, const std::complex<T> &value)
{
	std::valarray<std::complex<T> > Result(varray.size());
	for (int n=0; n<varray.size();n++)
		Result[n]=varray[n]+value;

	return Result;
}

template <class T>
inline std::valarray<std::complex<T> > operator+(const std::complex<T> &value, const std::valarray<T> &varray) {return varray+value;}

//Multiplication operator
template <class T>
inline std::valarray<std::complex<T> > operator*(const std::valarray<std::complex<T> > &varray1, const std::valarray<T> &varray2)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
		if (varray1.size()!=varray2.size()) throw("operator* valarray operands must have the same size.");
#endif

	std::valarray<std::complex<T> > Result(varray1);
	for (int n=0; n<varray1.size();n++)
		Result[n]*=varray2[n];

	return Result;
}

template <class T>
inline std::valarray<std::complex<T> > operator*(const std::valarray<T> &varray2, const std::valarray<std::complex<T> > &varray1) {return varray1*varray2;}

template <class T>
inline std::valarray<std::complex<T> > operator*(const std::valarray<std::complex<T> > &varray, const T &value)
{
	std::valarray<std::complex<T> > Result(varray);
	Result*=std::complex<T>(value,0);

	return Result;
}

template <class T>
inline std::valarray<std::complex<T> > operator*(const T &value, const std::valarray<std::complex<T> > &varray) {return varray*value;}

template <class T>
inline std::valarray<std::complex<T> > operator*(const std::valarray<T> &varray, const std::complex<T> &value)
{
	std::valarray<std::complex<T> > Result(varray.size());
	for (int n=0; n<varray.size();n++)
		Result[n]=varray[n]*value;

	return Result;
}

template <class T>
inline std::valarray<std::complex<T> > operator*(const std::complex<T> &value, const std::valarray<T> &varray) {return varray*value;}

//Subtraction operator
template <class T>
inline std::valarray<std::complex<T> > operator-(const std::valarray<std::complex<T> > &varray1, const std::valarray<T> &varray2)
{
	return varray1+(-varray2);
}

template <class T>
inline std::valarray<std::complex<T> > operator-(const std::valarray<T> &varray1, const std::valarray<std::complex<T> > &varray2)
{
	return -varray2+varray1;
}

template <class T>
inline std::valarray<std::complex<T> > operator-(const std::valarray<std::complex<T> > &varray, const T &value)
{
	return varray+(-value);
}

template <class T>
inline std::valarray<std::complex<T> > operator-(const T &value, const std::valarray<std::complex<T> > &varray)
{
	return value+(-varray);
}

template <class T>
inline std::valarray<std::complex<T> > operator-(const std::valarray<T> &varray, const std::complex<T> &value)
{
	return varray+(-value);
}

template <class T>
inline std::valarray<std::complex<T> > operator-(const std::complex<T> &value, const std::valarray<T> &varray)
{
	return value+(-varray);
}

//Division operator
template <class T>
inline std::valarray<std::complex<T> > operator/(const std::valarray<std::complex<T> > &varray1, const std::valarray<T> &varray2)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
		if (varray1.size()!=varray2.size()) throw("operator/ valarray operands must have the same size.");
#endif

	std::valarray<std::complex<T> > Result(varray1.size());
	for (int n=0; n<varray1.size();n++)
		Result[n]=varray1[n]/varray2[n];

	return Result;
}

template <class T>
inline std::valarray<std::complex<T> > operator/(const std::valarray<T> &varray1, const std::valarray<std::complex<T> > &varray2)
{
#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
		if (varray1.size()!=varray2.size()) throw("operator/ valarray operands must have the same size.");
#endif

	std::valarray<std::complex<T> > Result(varray1.size());
	for (int n=0; n<varray1.size();n++)
		Result[n]=varray1[n]/varray2[n];

	return Result;
}

template <class T>
inline std::valarray<std::complex<T> > operator/(const std::valarray<std::complex<T> > &varray, const T &value)
{
	std::valarray<std::complex<T> > Result(varray);
	Result/=std::complex<T>(value,0);

	return Result;
}

template <class T>
inline std::valarray<std::complex<T> > operator/(const T &value, const std::valarray<std::complex<T> > &varray)
{
	std::valarray<std::complex<T> > Result(complex<T>(value,0),varray.size());
	Result/=varray;

	return Result;
}

template <class T>
inline std::valarray<std::complex<T> > operator/(const std::valarray<T> &varray, const std::complex<T> &value)
{
	std::valarray<std::complex<T> > Result(varray.size());
	for (int n=0; n<varray.size();n++)
		Result[n]=varray[n]/value;

	return Result;
}

template <class T>
inline std::valarray<std::complex<T> > operator/(const std::complex<T> &value, const std::valarray<T> &varray)
{
	std::valarray<std::complex<T> > Result(varray.size());
	for (int n=0; n<varray.size();n++)
		Result[n]=value/varray[n];

	return Result;
}

namespace std{

template <class T>
std::valarray<T> real(const std::valarray<std::complex<T> >& complexVal)
{
	size_t N=complexVal.size();
	std::valarray<T> realVec(N);
	for (size_t i=0; i<N; i++)
		realVec[i]=complexVal[i].real();

	return realVec;
}

template <class T>
std::valarray<T> imag(const std::valarray<std::complex<T> >& complexVal)
{
	size_t N=complexVal.size();
	std::valarray<T> imagVec(N);
	for (size_t i=0; i<N; i++)
		imagVec[i]=complexVal[i].imag();

	return imagVec;
}

template <class T>
std::valarray<std::complex<T> > conj(const std::valarray<std::complex<T> >& complexVal)
{
	size_t N=complexVal.size();
	std::valarray<std::complex<T> > conjVec(N);
	for (size_t i=0; i<N; i++)
		conjVec[i]=std::conj(complexVal[i]);

	return conjVec;
}

template <class T>
std::valarray<T> abs(const std::valarray<std::complex<T> >& complexVal)
{
	size_t N=complexVal.size();
	std::valarray<T> absVec(N);
	for (size_t i=0; i<N; i++)
		absVec[i]=std::abs(complexVal[i]);

	return absVec;
}

template <class T>
std::valarray<T> arg(const std::valarray<std::complex<T> >& complexVal)
{
	size_t N=complexVal.size();
	std::valarray<T> argVec(N);
	for (size_t i=0; i<N; i++)
		argVec[i]=std::arg(complexVal[i]);

	return argVec;
}

}//end of the namespace std

//To overcome a bug in VS compiler
template <class T>
inline std::valarray<std::complex<T> > sqrt(const std::valarray<std::complex<T> >& ValArray)
{
	size_t N=ValArray.size();
	std::valarray<std::complex<T> > Res(N);
	for (size_t n=0;n<N;n++)
		Res[n]=std::sqrt(ValArray[n]);

	return Res;
}

#endif