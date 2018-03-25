#include <fintrf.h>

MODULE MX_INTERFACES
INTERFACE
	mwPointer function mxGetM(pm)
		mwPointer pm
	end function mxGetM

	mwPointer function mxGetN(pm)
		mwPointer pm
	end function mxGetN

	mwPointer function mxGetPr(pm)
		mwPointer pm
	end function mxGetPr

	mwPointer function mxGetPi(pm)
		mwPointer pm
	end function mxGetPi

	mwPointer function mxCreateDoubleMatrix(m, n, ComplexFlag)
		mwPointer m, n
        integer*4 ComplexFlag
	end function mxCreateDoubleMatrix

	mwPointer function mxCreateNumericMatrix(m, n, classid, ComplexFlag)
		mwSize m, n
        integer*4 classid, ComplexFlag
	end function mxCreateNumericMatrix

    real*8 function mxGetScalar(pm)
        mwPointer pm
    end function mxGetScalar

    subroutine mxCopyPtrToReal8(px, y, n)
        mwPointer px
        real*8 y(n)
        mwSize n
    end subroutine mxCopyPtrToReal8

	subroutine mxCopyReal8ToPtr(y, px, n)
        real*8 y(n)
        mwPointer px
        mwSize n
	end subroutine mxCopyReal8ToPtr

	subroutine mxCopyReal4ToPtr(y, px, n)
        real*4 y(n)
        mwPointer px
        mwSize n
	end subroutine mxCopyReal4ToPtr

	subroutine mxCopyPtrToComplex16(pr, pi, y, n)
        mwPointer pr, pi
        complex*16 y(n)
        mwSize n
	end subroutine mxCopyPtrToComplex16

	subroutine mxCopyPtrToComplex8(pr, pi, y, n)
        mwPointer pr, pi
        complex*8 y(n)
        mwSize n
	end subroutine mxCopyPtrToComplex8

	subroutine mxCopyComplex16ToPtr(y, pr, pi, n)
        complex*16 y(n)
        mwPointer pr, pi
        mwSize n
	end subroutine mxCopyComplex16ToPtr

	subroutine mxCopyComplex8ToPtr(y, pr, pi, n)
        complex*8 y(n)
        mwPointer pr, pi
        mwSize n
	end subroutine mxCopyComplex8ToPtr

	integer*4 function mxIsComplex(pm)
		mwPointer pm
	end function mxIsComplex

	integer*4 function mxIsSingle(pm)
		mwPointer pm
	end function mxIsSingle

	integer*4 function mxIsDouble(pm)
		mwPointer pm
	end function mxIsDouble

	integer*4 function mxClassIDFromClassName(classname)
		character*(*) classname
	end function mxClassIDFromClassName

	mwPointer function mxGetNumberOfElements(pm)
		mwPointer pm
	end function mxGetNumberOfElements

	integer*4 function mxGetClassID(pm)
		mwPointer pm
	end function mxGetClassID

    integer*4 function mxIsLogical(pm)
        mwPointer pm
    end function mxIsLogical

    subroutine mxCopyPtrToInteger1(px, y, n)
        mwPointer px
        integer*1 y(n)
        mwSize n
    end subroutine

    mwPointer function mxGetData(pm)
        mwPointer pm
    end function

END INTERFACE
END MODULE MX_INTERFACES