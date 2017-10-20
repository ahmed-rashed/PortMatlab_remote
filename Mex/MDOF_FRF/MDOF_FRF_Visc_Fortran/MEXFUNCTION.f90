#include <fintrf.h>

SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
! This subroutine is the main gateway to MATLAB.  When a MEX function
!  is executed MATLAB calls this subroutine

USE IFPORT
USE MX_INTERFACES
USE MEX_INTERFACES

IMPLICIT NONE
!DEC$ ATTRIBUTES DLLEXPORT :: MEXFUNCTION
MWPOINTER PLHS(*), PRHS(*)
INTEGER(4) NLHS, NRHS
!mwSize NN
mwPointer N, N_w, N_cols

INTEGER(4), ALLOCATABLE :: n_vec(:), m_vec(:)
REAL(8), ALLOCATABLE :: n_vec_temp(:), m_vec_temp(:)
COMPLEX(8), ALLOCATABLE :: EigValues_vec(:), EigVectors_Normalized_mat(:,:), H_cols(:,:)
REAL(8), ALLOCATABLE :: w_column(:,:)
INTEGER(4) ii

INTERFACE
    FUNCTION MDOF_FRF_Visc(EigValues_vec, EigVectors_Normalized_mat, N, w_column, N_w, n_vec, m_vec, N_cols)
        INTEGER(4), INTENT(IN) :: N, N_w, N_cols
        INTEGER(4), INTENT(IN) :: n_vec(N_cols), m_vec(N_cols)
        COMPLEX(8), INTENT(IN) :: EigValues_vec(2*N), EigVectors_Normalized_mat(N,2*N)
        REAL(8), INTENT(IN) :: w_column(N_w,1)
        COMPLEX(8) MDOF_FRF_Visc(N_w,N_cols)
    END FUNCTION MDOF_FRF_Visc
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
IF ((NRHS/=5)) THEN
	CALL mexErrMsgTxt('Dear student, This function needs 5 inputs.')
END IF

IF ((NLHS>1)) THEN
	CALL mexErrMsgTxt('Dear student, this function returns only one output.')
END IF

N=mxGetM(PRHS(2))
IF (2*N /= mxGetN(PRHS(2))) THEN
	CALL mexErrMsgTxt('Dear student, The EigVectors_Normalized matrix must have columns twice as many as be rows.')
END IF

IF (mxGetNumberOfElements(prhs(1)) /= 2*N) THEN
	CALL mexErrMsgTxt('Dear student, the EigValues_vec vector must have length identical to the number of columns of EigVectors_Normalized matrix.')
END IF

IF (mxGetN(prhs(3)) /= 1) THEN
	CALL mexErrMsgTxt('Dear student, w_column must be a column vector.')
END IF
N_w=mxGetM(prhs(3))

IF (mxGetM(prhs(4)) /= 1) THEN
	CALL mexErrMsgTxt('Dear student, n_row must be a row vector.')
END IF
N_cols=mxGetN(prhs(4))

IF ((mxGetM(prhs(5)) /= 1) .OR. (mxGetN(prhs(5)) /= N_cols)) THEN
	CALL mexErrMsgTxt('Dear student, m_row must have size identical to n_row.')
END IF

!Inputs
ALLOCATE(EigValues_vec(2*N), EigVectors_Normalized_mat(N,2*N),w_column(N_w,1), n_vec(N_cols), m_vec(N_cols), n_vec_temp(N_cols), m_vec_temp(N_cols))
CALL mxCopyPtrToComplex16(mxGetPr(PRHS(1)), mxGetPi(PRHS(1)), EigValues_vec, 2*N)
CALL mxCopyPtrToComplex16(mxGetPr(PRHS(2)), mxGetPi(PRHS(2)), EigVectors_Normalized_mat, 2*N*N)
CALL mxCopyPtrToReal8(mxGetPr(PRHS(3)), w_column, N_w)
CALL mxCopyPtrToReal8(mxGetPr(PRHS(4)), n_vec_temp, N_cols)
CALL mxCopyPtrToReal8(mxGetPr(PRHS(5)), m_vec_temp, N_cols)

DO ii=1,N_cols
    n_vec=n_vec_temp
    m_vec=m_vec_temp
END DO
!Do Actual Computations
ALLOCATE(H_cols(N_w,N_cols))
H_cols=MDOF_FRF_Visc(EigValues_vec, EigVectors_Normalized_mat, N, w_column, N_w, n_vec, m_vec, N_cols)

PLHS(1)=mxCreateDoubleMatrix(N_w,N_cols,1)  !Create the output matrix
CALL mxCopyComplex16ToPtr(H_cols, mxGetPr(PLHS(1)), mxGetPi(PLHS(1)), N_w*N_cols)

END SUBROUTINE mexFunction