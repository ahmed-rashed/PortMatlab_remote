#include <fintrf.h>

logical function mxGetLogicalScalar(mx)
    USE MX_INTERFACES

    implicit none
    mwPointer, intent(in) :: mx

    integer*1 integer1value(1)

    mxGetLogicalScalar = .false.
    IF ((mxIsLogical(mx) /= 0) .AND. (mxGetNumberOfElements(mx) == 1)) THEN
        CALL mxCopyPtrToInteger1(mxGetData(mx),integer1value,1)
        mxGetLogicalScalar= (integer1value(1) /= 0)
    END IF
end function

SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
! This subroutine is the main gateway to MATLAB.  When a MEX function
!  is executed MATLAB calls this subroutine

USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_SIZE_T
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
    !mwSize NN
    mwPointer N

    REAL(8), ALLOCATABLE :: M_mat(:,:),C_mat(:,:),K_mat(:,:)
    LOGICAL isProportional
    INTEGER(C_SIZE_T) NN
    COMPLEX(8), ALLOCATABLE :: EigVectors_Normalized(:,:), EigValues_vec(:)

    INTERFACE
        SUBROUTINE MDOF_Eig_Visc(M_mat,C_mat,K_mat,N,isProportional,EigValues_vec,EigVectors_Normalized)
            USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_SIZE_T
            INTEGER(C_SIZE_T), INTENT(IN) :: N
            LOGICAL, INTENT(IN) :: isProportional
            REAL(8), INTENT(IN) ::  M_mat(N,N),C_mat(N,N),K_mat(N,N)
            COMPLEX(8), INTENT(OUT) :: EigVectors_Normalized(N,2*N), EigValues_vec(2*N)
	    END SUBROUTINE MDOF_Eig_Visc

        logical function mxGetLogicalScalar(mx)
            mwPointer, intent(in) :: mx
        end function mxGetLogicalScalar
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
    IF ((NRHS/=3) .AND. (NRHS/=4)) THEN
	    CALL mexErrMsgTxt('Dear student, This function takes 3 or 4 inputs.')
    END IF

    IF ((NLHS/=1) .AND. (NLHS/=2)) THEN
	    CALL mexErrMsgTxt('Dear student, this function returns one or two outputs only.')
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

    IF (NRHS==4) THEN
        IF ((mxGetNumberOfElements(PRHS(4)) /= 1) .OR. (.NOT. mxIsLogical(PRHS(4)))) THEN
            CALL mexErrMsgTxt('Dear student, isProportional must be a logical scalar/integer.')
        END IF
    END IF

    !Inputs
    ALLOCATE(M_mat(N,N),C_mat(N,N),K_mat(N,N))
    CALL mxCopyPtrToReal8(mxGetDoubles(PRHS(1)), M_mat, N*N)
    CALL mxCopyPtrToReal8(mxGetDoubles(PRHS(2)), C_mat, N*N)
    CALL mxCopyPtrToReal8(mxGetDoubles(PRHS(3)), K_mat, N*N)
    isProportional=.FALSE.
    IF (NRHS==4) THEN
        isProportional=mxGetLogicalScalar(PRHS(4))
    END IF


    !Do Actual Computations
    ALLOCATE(EigVectors_Normalized(N,2*N), EigValues_vec(2*N))
    CALL MDOF_Eig_Visc(M_mat,C_mat,K_mat,N,isProportional,EigValues_vec,EigVectors_Normalized)

    !Outputs
    PLHS(1)=mxCreateDoubleMatrix(N,2*N,1)  !Create the output matrix
    CALL mxCopyComplex16ToPtr(EigVectors_Normalized, mxGetComplexDoubles(PLHS(1)), N*2*N)

    IF (NLHS > 1) THEN
        PLHS(2)=mxCreateDoubleMatrix(2*N,1,1)  ! Create the output matrix
    CALL mxCopyComplex16ToPtr(EigValues_vec, mxGetComplexDoubles(PLHS(2)), 2*N)
END IF

END SUBROUTINE MEXFUNCTION