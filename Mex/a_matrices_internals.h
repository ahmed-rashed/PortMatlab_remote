#ifndef a_matrices_internals_h
#define a_matrices_internals_h

#include <iostream>

namespace a_matrices{

#if (defined(_CHECKBOUNDS_) || defined(_DEBUG))
template <class T>
class next_sub
{
public:
	const T* prow;
	const size_t cols;

	//constructor
	next_sub(const T* pprow,const size_t& ccols):prow(pprow), cols(ccols) {}

	T& operator[](const size_t &j)
	{
		if (j>=cols) throw("Column subscript is out of bounds.");
		return (const_cast<T*>(prow))[j];
	}

	const T& operator[](const size_t &j) const
	{
		if (j>=cols) throw("Column subscript is out of bounds.");
		return (const_cast<T*>(prow))[j];
	}
};
#endif


}//end of the namespace a_matrices
#endif