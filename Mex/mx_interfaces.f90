#include <fintrf.h>

MODULE MX_INTERFACES
INTERFACE
	mwPointer function mxGetM(pm)
	    IMPLICIT NONE
		mwPointer pm
	end function mxGetM

	mwPointer function mxGetN(pm)
	    IMPLICIT NONE
		mwPointer pm
	end function mxGetN

	mwPointer function mxGetDoubles(pm)
		IMPLICIT NONE
		mwPointer pm
	end function mxGetDoubles

	mwPointer function mxGetComplexDoubles(pm)
		IMPLICIT NONE
		mwPointer pm
	end function mxGetComplexDoubles

	mwPointer function mxCreateDoubleMatrix(m, n, ComplexFlag)
		IMPLICIT NONE
		mwPointer m, n
        integer*4 ComplexFlag
	end function mxCreateDoubleMatrix

	mwPointer function mxCreateNumericMatrix(m, n, classid, ComplexFlag)
		IMPLICIT NONE
		mwSize m, n
        integer*4 classid, ComplexFlag
	end function mxCreateNumericMatrix

    real*8 function mxGetScalar(pm)
		IMPLICIT NONE
        mwPointer pm
    end function mxGetScalar

    subroutine mxCopyPtrToReal8(px, y, n)
        IMPLICIT NONE
        mwPointer px
        real*8 y(n)
        mwSize n
    end subroutine mxCopyPtrToReal8

	subroutine mxCopyReal8ToPtr(y, px, n)
	    IMPLICIT NONE
        real*8 y(n)
        mwPointer px
        mwSize n
	end subroutine mxCopyReal8ToPtr

	subroutine mxCopyReal4ToPtr(y, px, n)
	    IMPLICIT NONE
        real*4 y(n)
        mwPointer px
        mwSize n
	end subroutine mxCopyReal4ToPtr

	subroutine mxCopyPtrToComplex16(pd, y, n)
	    IMPLICIT NONE
        mwPointer pd
        complex*16 y(n)
        mwSize n
	end subroutine mxCopyPtrToComplex16

	subroutine mxCopyPtrToComplex8(pd, y, n)
	    IMPLICIT NONE
        mwPointer pd
        complex*8 y(n)
        mwSize n
	end subroutine mxCopyPtrToComplex8

	subroutine mxCopyComplex16ToPtr(y, pd, n)
	    IMPLICIT NONE
        complex*16 y(n)
        mwPointer pd
        mwSize n
	end subroutine mxCopyComplex16ToPtr

	subroutine mxCopyComplex8ToPtr(y, pd, n)
	    IMPLICIT NONE
        complex*8 y(n)
        mwPointer pd
        mwSize n
	end subroutine mxCopyComplex8ToPtr

	integer*4 function mxIsComplex(pm)
		IMPLICIT NONE
		mwPointer pm
	end function mxIsComplex

	integer*4 function mxIsSingle(pm)
		IMPLICIT NONE
		mwPointer pm
	end function mxIsSingle

	integer*4 function mxIsDouble(pm)
		IMPLICIT NONE
		mwPointer pm
	end function mxIsDouble

	integer*4 function mxClassIDFromClassName(classname)
		IMPLICIT NONE
		character*(*) classname
	end function mxClassIDFromClassName

	mwPointer function mxGetNumberOfElements(pm)
		IMPLICIT NONE
		mwPointer pm
	end function mxGetNumberOfElements

	integer*4 function mxGetClassID(pm)
		IMPLICIT NONE
		mwPointer pm
	end function mxGetClassID

    integer*4 function mxIsLogical(pm)
		IMPLICIT NONE
        mwPointer pm
    end function mxIsLogical

    subroutine mxCopyPtrToInteger1(px, y, n)
		IMPLICIT NONE
        mwPointer px
        integer*1 y(n)
        mwSize n
    end subroutine

    mwPointer function mxGetData(pm)
		IMPLICIT NONE
        mwPointer pm
    end function

END INTERFACE
END MODULE MX_INTERFACES
