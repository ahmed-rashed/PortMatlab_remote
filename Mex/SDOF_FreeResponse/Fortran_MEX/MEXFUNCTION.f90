#include <fintrf.h>

SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS) bind(C, name="MEXFUNCTION")
! This subroutine is the main gateway to MATLAB.  When a MEX function
!  is executed MATLAB calls this subroutine

!USE IFPORT
USE MX_INTERFACES
USE MEX_INTERFACES
USE iso_c_binding, only: c_int
use iso_c_binding, only: C_CHAR, C_NULL_CHAR

IMPLICIT NONE
!DEC$ ATTRIBUTES DLLEXPORT :: MEXFUNCTION
MWPOINTER PLHS(*), PRHS(*)
INTEGER(c_int), VALUE :: NLHS, NRHS
REAL(8) w_n,zeta,x0,v0
mwSize N
REAL(8), ALLOCATABLE :: t(:),x(:)

INTERFACE
    FUNCTION FREE_RESPONSE(w_n,zeta,x0,v0,t,N)
        INTEGER(4) N
        REAL(8) w_n,zeta,x0,v0,t(N),FREE_RESPONSE(N)
    END FUNCTION FREE_RESPONSE
END INTERFACE

#if defined(_WIN32) || defined(_WIN64)
! For Windows only!
! This resets the floating point exception to allow divide by zero,
! overflow and invalid numbers.
!
    INTEGER(2) CONTROL
    CALL GETCONTROLFPQQ(CONTROL)
    CONTROL = CONTROL .OR. FPCW$ZERODIVIDE
      CONTROL = CONTROL .OR. FPCW$INVALID
      CONTROL = CONTROL .OR. FPCW$OVERFLOW
    CALL SETCONTROLFPQQ(CONTROL)
#endif

! Check for proper number of arguments
IF ((NRHS/=5)) THEN
    CALL mexErrMsgTxt(C_CHAR_"Dear student, This function needs 5 inputs."//C_NULL_CHAR)
ENDIF

IF ((NLHS>1)) THEN
    CALL mexErrMsgTxt(C_CHAR_"Dear student, this function returns only one output."//C_NULL_CHAR)
ENDIF

DO N=1,4
    IF (mxGetNumberOfElements(PRHS(N)) /= 1) THEN
        CALL mexErrMsgTxt(C_CHAR_"Dear student, first four inputs must be scalars."//C_NULL_CHAR)    !Improve this
    ENDIF
END DO

w_n=mxGetScalar(PRHS(1))
zeta=mxGetScalar(PRHS(2))
x0=mxGetScalar(PRHS(3))
v0=mxGetScalar(PRHS(4))
N=mxGetNumberOfElements(PRHS(5))
ALLOCATE(t(N),x(N))
CALL mxCopyPtrToReal8(mxGetPr(PRHS(5)), t, N)

x=FREE_RESPONSE(w_n,zeta,x0,v0,t,N)

PLHS(1)=mxCreateDoubleMatrix(mxGetM(PRHS(5)),mxGetN(PRHS(5)),0)
CALL mxCopyReal8ToPtr(x, mxGetPr(PLHS(1)), N)

END SUBROUTINE MEXFUNCTION
