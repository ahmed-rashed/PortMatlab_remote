#include <fintrf.h>

SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
! This subroutine is the main gateway to MATLAB.  When a MEX function
!  is executed MATLAB calls this subroutine

USE IFPORT
USE MX_INTERFACES
USE MEX_INTERFACES

IMPLICIT NONE
#if defined(_WIN32) || defined(_WIN64)
!DEC$ ATTRIBUTES DLLEXPORT :: MEXFUNCTION
#else
!GCC$ ATTRIBUTES DLLEXPORT :: MEXFUNCTION
#endif
MWPOINTER PLHS(*), PRHS(*)
INTEGER(4) NLHS, NRHS
mwSize NN
mwPointer N,n_w_points,n_n_row

REAL(8), ALLOCATABLE :: M_mat(:,:),C_mat(:,:),K_mat(:,:),w_col(:,:),n_row(:,:),m_row(:,:)
COMPLEX(8), ALLOCATABLE :: H_cols(:,:)

INTERFACE
    FUNCTION MDOF_FRF_VISC_SLOW(M_mat,C_mat,K_mat,N,w_col,n_w_points,n_vec,m_vec,n_n_row)
        USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_SIZE_T
        INTEGER(C_SIZE_T), INTENT(IN) :: N,n_w_points,n_n_row
        REAL(8), INTENT(IN) ::  M_mat(N,N),C_mat(N,N),K_mat(N,N),w_col(n_w_points,1),n_vec(n_n_row),m_vec(n_n_row)
        COMPLEX(8) MDOF_FRF_VISC_SLOW(n_w_points,n_n_row)
    END FUNCTION MDOF_FRF_VISC_SLOW
END INTERFACE

#if defined(_WIN32) || defined(_WIN64)
! For Windows only!
! This resets the floating point exception to allow divide by zero,
! overflow and invalid numbers. 
    INTEGER(2) CONTROL
    CALL GETCONTROLFPQQ(CONTROL)
    CONTROL = CONTROL .OR. FPCW$ZERODIVIDE
      CONTROL = CONTROL .OR. FPCW$INVALID
      CONTROL = CONTROL .OR. FPCW$OVERFLOW
    CALL SETCONTROLFPQQ(CONTROL)
#endif

! Check for proper number of arguments
IF ((NRHS/=6)) THEN
    CALL mexErrMsgTxt('Dear student, This function needs 6 inputs.')
END IF

IF ((NLHS>1)) THEN
    CALL mexErrMsgTxt('Dear student, this function returns only one output.')
END IF

N=mxGetM(PRHS(1))
IF (N /= mxGetN(PRHS(1))) THEN
    CALL mexErrMsgTxt('Dear student, The M_mat matrix must be square.')
END IF

DO NN=2,3
    IF ((mxGetM(PRHS(NN)) /= N) .OR. (mxGetN(PRHS(NN)) /= N)) THEN
        CALL mexErrMsgTxt('Dear student, the C_mat and K_mat matrices must be square with size identical to M_mat.')    !Improve this
    END IF
END DO

IF (mxGetN(PRHS(4)) /= 1) THEN
    CALL mexErrMsgTxt('Dear student, w_col must be a column vector.')
END IF
n_w_points=mxGetM(PRHS(4));

IF (mxGetM(PRHS(5)) /= 1) THEN
    CALL mexErrMsgTxt('Dear student, n_row must be a row vector.')
END IF
n_n_row=mxGetN(PRHS(5))

IF (mxGetM(PRHS(6)) /= 1) THEN
    CALL mexErrMsgTxt('Dear student, m_row must be a row vector.')
END IF

IF (mxGetN(PRHS(6)) /= n_n_row) THEN
    CALL mexErrMsgTxt('Dear student, m_row & m_row must have the same size.')
END IF

!Inputs
ALLOCATE(M_mat(N,N),C_mat(N,N),K_mat(N,N),w_col(n_w_points,1),n_row(1,n_n_row),m_row(1,n_n_row))
CALL mxCopyPtrToReal8(mxGetDoubles(PRHS(1)), M_mat, N*N)
CALL mxCopyPtrToReal8(mxGetDoubles(PRHS(2)), C_mat, N*N)
CALL mxCopyPtrToReal8(mxGetDoubles(PRHS(3)), K_mat, N*N)
CALL mxCopyPtrToReal8(mxGetDoubles(PRHS(4)), w_col, n_w_points)
CALL mxCopyPtrToReal8(mxGetDoubles(PRHS(5)), n_row, n_n_row)
CALL mxCopyPtrToReal8(mxGetDoubles(PRHS(6)), m_row, n_n_row)

!Do Actual Computations
ALLOCATE(H_cols(n_w_points,n_n_row))
H_cols=MDOF_FRF_VISC_SLOW(M_mat,C_mat,K_mat,N,w_col,n_w_points,n_row,m_row,n_n_row)

!Outputs
PLHS(1)=mxCreateDoubleMatrix(n_w_points,n_n_row,1)  ! Create the output matrix
CALL mxCopyComplex16ToPtr(H_cols, mxGetComplexDoubles(PLHS(1)), n_w_points*n_n_row)

END SUBROUTINE MEXFUNCTION
