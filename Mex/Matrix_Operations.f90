Function sMatDiag(A_mat,N)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_SIZE_T
    
	INTEGER(C_SIZE_T), INTENT(IN) :: N
	REAL(4), INTENT(IN) :: A_mat(N,N)
	REAL(4) sMatDiag(N)

    INTEGER(C_SIZE_T) NN

    DO NN=1,N
        sMatDiag(NN)=A_mat(NN,NN)
    END DO
END FUNCTION sMatDiag

Function dMatDiag(A_mat,N)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_SIZE_T
    
	INTEGER(C_SIZE_T), INTENT(IN) :: N
	REAL(8), INTENT(IN) :: A_mat(N,N)
	REAL(8) dMatDiag(N)

    INTEGER(C_SIZE_T) NN

    DO NN=1,N
        dMatDiag(NN)=A_mat(NN,NN)
    END DO
END FUNCTION dMatDiag

Function cMatDiag(A_mat,N)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_SIZE_T
    
	INTEGER(C_SIZE_T), INTENT(IN) :: N
	COMPLEX(4), INTENT(IN) :: A_mat(N,N)
	COMPLEX(4) cMatDiag(N)

    INTEGER(C_SIZE_T) NN

    DO NN=1,N
        cMatDiag(NN)=A_mat(NN,NN)
    END DO
END FUNCTION cMatDiag

Function zMatDiag(A_mat,N)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_SIZE_T
    
	INTEGER(C_SIZE_T), INTENT(IN) :: N
	COMPLEX(8), INTENT(IN) :: A_mat(N,N)
	COMPLEX(8) zMatDiag(N)

    INTEGER(C_SIZE_T) NN

    DO NN=1,N
        zMatDiag(NN)=A_mat(NN,NN)
    END DO
END FUNCTION zMatDiag

Function sDiagMat(VEC,N)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_SIZE_T
    
	INTEGER(C_SIZE_T), INTENT(IN) :: N
	REAL(4), INTENT(IN) :: VEC(N)
	REAL(4) sDiagMat(N,N)

    INTEGER(C_SIZE_T) NN

    sDiagMat=0.0
    DO NN=1,N
        sDiagMat(NN,NN)=VEC(NN)
    END DO
END FUNCTION sDiagMat

Function dDiagMat(VEC,N)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_SIZE_T
    
	INTEGER(C_SIZE_T), INTENT(IN) :: N
	REAL(8), INTENT(IN) :: VEC(N)
	REAL(8) dDiagMat(N,N)

    INTEGER(C_SIZE_T) NN

    dDiagMat=0.0
    DO NN=1,N
        dDiagMat(NN,NN)=VEC(NN)
    END DO
END FUNCTION dDiagMat

Function cDiagMat(VEC,N)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_SIZE_T
    
	INTEGER(C_SIZE_T), INTENT(IN) :: N
	COMPLEX(4), INTENT(IN) :: VEC(N)
	COMPLEX(4) cDiagMat(N,N)

    INTEGER(C_SIZE_T) NN

    cDiagMat=0.0
    DO NN=1,N
        cDiagMat(NN,NN)=VEC(NN)
    END DO
END FUNCTION cDiagMat

Function zDiagMat(VEC,N)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_SIZE_T
    
	INTEGER(C_SIZE_T), INTENT(IN) :: N
	COMPLEX(8), INTENT(IN) :: VEC(N)
	COMPLEX(8) zDiagMat(N,N)

    INTEGER(C_SIZE_T) NN

    zDiagMat=0.0
    DO NN=1,N
        zDiagMat(NN,NN)=VEC(NN)
    END DO
END FUNCTION zDiagMat