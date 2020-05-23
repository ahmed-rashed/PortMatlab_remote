#include <fintrf.h>

SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
    ! This subroutine is the main gateway to MATLAB.  When a MEX function
    !  is executed MATLAB calls this subroutine

    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_SIZE_T
#if defined(_WIN32) || defined(_WIN64)
    USE IFPORT
#endif
    USE MX_INTERFACES
    USE MEX_INTERFACES

    IMPLICIT NONE
#if defined(_WIN32) || defined(_WIN64)
!DEC$ ATTRIBUTES DLLEXPORT :: MEXFUNCTION
#endif
    MWPOINTER PLHS(*), PRHS(*)
    INTEGER(4) NLHS, NRHS
    REAL(8) w_n,zeta,x0,v0
    mwSize N
    INTEGER(C_SIZE_T) NN
    REAL(8), ALLOCATABLE :: t(:),x(:)

    INTERFACE
        FUNCTION FREE_RESPONSE(w_n,zeta,x0,v0,t,N)
            USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_SIZE_T
            INTEGER(C_SIZE_T), INTENT(IN) :: N
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
        CALL mexErrMsgTxt('Dear student, This function needs 5 inputs.')
    ENDIF

    IF ((NLHS>1)) THEN
        CALL mexErrMsgTxt('Dear student, this function returns only one output.')
    ENDIF

    DO N=1,4
        IF (mxGetNumberOfElements(PRHS(N)) /= 1) THEN
            CALL mexErrMsgTxt('Dear student, first four inputs must be scalars.')    !Improve this
        ENDIF
    END DO

    w_n=mxGetScalar(PRHS(1))
    zeta=mxGetScalar(PRHS(2))
    x0=mxGetScalar(PRHS(3))
    v0=mxGetScalar(PRHS(4))
    N=mxGetNumberOfElements(PRHS(5))
    ALLOCATE(t(N),x(N))
    CALL mxCopyPtrToReal8(mxGetDoubles(PRHS(5)), t, N)

    NN=N
    x=FREE_RESPONSE(w_n,zeta,x0,v0,t,NN)

    PLHS(1)=mxCreateDoubleMatrix(mxGetM(PRHS(5)),mxGetN(PRHS(5)),INT(0,4))
    CALL mxCopyReal8ToPtr(x, mxGetDoubles(PLHS(1)), N)

END SUBROUTINE MEXFUNCTION
