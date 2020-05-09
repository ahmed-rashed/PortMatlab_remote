MODULE Matrix_Operations_interfaces

INTERFACE

    Function sMatDiag(A_mat,N)
        USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_SIZE_T
        
	    INTEGER(C_SIZE_T), INTENT(IN) :: N
	    REAL(4), INTENT(IN) :: A_mat(N,N)
	    REAL(4) sMatDiag(N)
    END FUNCTION sMatDiag

    Function dMatDiag(A_mat,N)
        USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_SIZE_T
        
	    INTEGER(C_SIZE_T), INTENT(IN) :: N
	    REAL(8), INTENT(IN) :: A_mat(N,N)
	    REAL(8) dMatDiag(N)
    END FUNCTION dMatDiag

    Function cMatDiag(A_mat,N)
        USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_SIZE_T
        
	    INTEGER(C_SIZE_T), INTENT(IN) :: N
	    COMPLEX(4), INTENT(IN) :: A_mat(N,N)
	    COMPLEX(4) cMatDiag(N)
    END FUNCTION cMatDiag

    Function zMatDiag(A_mat,N)
        USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_SIZE_T
        
	    INTEGER(C_SIZE_T), INTENT(IN) :: N
	    COMPLEX(8), INTENT(IN) :: A_mat(N,N)
	    COMPLEX(8) zMatDiag(N)
    END FUNCTION zMatDiag

    Function sDiagMat(VEC,N)
        USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_SIZE_T
        
	    INTEGER(C_SIZE_T), INTENT(IN) :: N
	    REAL(4), INTENT(IN) :: VEC(N)
	    REAL(4) sDiagMat(N,N)
    END FUNCTION sDiagMat

    Function dDiagMat(VEC,N)
        USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_SIZE_T
        
	    INTEGER(C_SIZE_T), INTENT(IN) :: N
	    REAL(8), INTENT(IN) :: VEC(N)
	    REAL(8) dDiagMat(N,N)
    END FUNCTION dDiagMat

    Function cDiagMat(VEC,N)
        USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_SIZE_T
        
	    INTEGER(C_SIZE_T), INTENT(IN) :: N
	    COMPLEX(4), INTENT(IN) :: VEC(N)
	    COMPLEX(4) cDiagMat(N,N)
    END FUNCTION cDiagMat

    Function zDiagMat(VEC,N)
        USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_SIZE_T
        
	    INTEGER(C_SIZE_T), INTENT(IN) :: N
	    COMPLEX(8), INTENT(IN) :: VEC(N)
	    COMPLEX(8) zDiagMat(N,N)
    END FUNCTION zDiagMat
END INTERFACE

END MODULE Matrix_Operations_interfaces