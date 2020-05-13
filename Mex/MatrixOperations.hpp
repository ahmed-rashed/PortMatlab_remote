#ifndef MatrixOperations_hpp
#define MatrixOperations_hpp

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
	for (size_t n=0; n<ind.size(); n++)
	{
		RowMajorInd2Sub(n_cols, ind[n], row[n], col[n]);
	}
}

template <template<typename> class T>
inline void ColMajorInd2Sub(const size_t& n_rows, const T<size_t>& ind, T<size_t>& row, T<size_t>& col)
{
	RowMajorInd2Sub(n_rows, ind, col, row);
}


#endif