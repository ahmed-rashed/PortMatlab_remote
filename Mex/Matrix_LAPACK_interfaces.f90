MODULE Matrix_LAPACK_interfaces

INTERFACE

    FUNCTION sInv(A,N)
        INTEGER(4), INTENT(IN) :: N
        REAL(4), INTENT(IN) :: A(N,N)
        REAL(4) sInv(N,N)
    END FUNCTION sInv

    FUNCTION dInv(A,N)
        INTEGER(4), INTENT(IN) :: N
        REAL(8), INTENT(IN) :: A(N,N)
        REAL(8) dInv(N,N)
    END FUNCTION dInv

    FUNCTION cInv(A,N)
        INTEGER(4), INTENT(IN) :: N
        COMPLEX(4), INTENT(IN) :: A(N,N)
        COMPLEX(4) cInv(N,N)
    END FUNCTION cInv

    FUNCTION zInv(A,N)
        INTEGER(4), INTENT(IN) :: N
        COMPLEX(8), INTENT(IN) :: A(N,N)
        COMPLEX(8) zInv(N,N)
    END FUNCTION zInv

    FUNCTION sSymInv(A_mat, N)
	    INTEGER(4), INTENT(IN) :: N
	    REAL(4), INTENT(IN) :: A_mat(N,N)
    END FUNCTION sSymInv

    FUNCTION dSymInv(A_mat, N)
	    INTEGER(4), INTENT(IN) :: N
	    REAL(8), INTENT(IN) :: A_mat(N,N)
    END FUNCTION dSymInv

    FUNCTION cSymInv(A_mat, N)
	    INTEGER(4), INTENT(IN) :: N
	    COMPLEX(4), INTENT(IN) :: A_mat(N,N)
    END FUNCTION cSymInv

    FUNCTION zSymInv(A_mat, N)
	    INTEGER(4), INTENT(IN) :: N
	    COMPLEX(8), INTENT(IN) :: A_mat(N,N)
    END FUNCTION zSymInv

    FUNCTION sSymPosDefInv(A_mat, N)
	    INTEGER(4), INTENT(IN) :: N
	    REAL(4), INTENT(IN) :: A_mat(N,N)
    END FUNCTION sSymPosDefInv

    FUNCTION dSymPosDefInv(A_mat, N)
	    INTEGER(4), INTENT(IN) :: N
	    REAL(8), INTENT(IN) :: A_mat(N,N)
    END FUNCTION dSymPosDefInv

    FUNCTION cHermPosDefInv(A_mat, N)
	    INTEGER(4), INTENT(IN) :: N
	    COMPLEX(4), INTENT(IN) :: A_mat(N,N)
    END FUNCTION cHermPosDefInv

    FUNCTION zHermPosDefInv(A_mat, N)
	    INTEGER(4), INTENT(IN) :: N
	    COMPLEX(8), INTENT(IN) :: A_mat(N,N)
    END FUNCTION zHermPosDefInv


    FUNCTION sSolve(A_mat, B_mat, N, nCols)
	    INTEGER(4), INTENT(IN) :: N, nCols
	    REAL(4), INTENT(IN) :: A_mat(N,N), B_mat(N,nCols)
    END FUNCTION sSolve

    FUNCTION dSolve(A_mat, B_mat, N, nCols)
	    INTEGER(4), INTENT(IN) :: N, nCols
	    REAL(8), INTENT(IN) :: A_mat(N,N), B_mat(N,nCols)
    END FUNCTION dSolve

    FUNCTION cSolve(A_mat, B_mat, N, nCols)
	    INTEGER(4), INTENT(IN) :: N, nCols
	    COMPLEX(4), INTENT(IN) :: A_mat(N,N), B_mat(N,nCols)
    END FUNCTION cSolve

    FUNCTION zSolve(A_mat, B_mat, N, nCols)
	    INTEGER(4), INTENT(IN) :: N, nCols
	    COMPLEX(8), INTENT(IN) :: A_mat(N,N), B_mat(N,nCols)
    END FUNCTION zSolve

    FUNCTION sSymPosDefSolve(A_mat, B_mat, N, nCols)
	    INTEGER(4), INTENT(IN) :: N, nCols
	    REAL(4), INTENT(IN) :: A_mat(N,N), B_mat(N,nCols)
    END FUNCTION sSymPosDefSolve

    FUNCTION dSymPosDefSolve(A_mat, B_mat, N, nCols)
	    INTEGER(4), INTENT(IN) :: N, nCols
	    REAL(8), INTENT(IN) :: A_mat(N,N), B_mat(N,nCols)
    END FUNCTION dSymPosDefSolve

    FUNCTION cHermPosDefSolve(A_mat, B_mat, N, nCols)
	    INTEGER(4), INTENT(IN) :: N, nCols
	    COMPLEX(4), INTENT(IN) :: A_mat(N,N), B_mat(N,nCols)
    END FUNCTION cHermPosDefSolve

    FUNCTION zHermPosDefSolve(A_mat, B_mat, N, nCols)
	    INTEGER(4), INTENT(IN) :: N, nCols
	    COMPLEX(8), INTENT(IN) :: A_mat(N,N), B_mat(N,nCols)
    END FUNCTION zHermPosDefSolve

    FUNCTION sSymSolve(A_mat, B_mat, N, nCols)
	    INTEGER(4), INTENT(IN) :: N, nCols
	    REAL(4), INTENT(IN) :: A_mat(N,N), B_mat(N,nCols)
    END FUNCTION sSymSolve

    FUNCTION dSymSolve(A_mat, B_mat, N, nCols)
	    INTEGER(4), INTENT(IN) :: N, nCols
	    REAL(8), INTENT(IN) :: A_mat(N,N), B_mat(N,nCols)
    END FUNCTION dSymSolve

    FUNCTION cSymSolve(A_mat, B_mat, N, nCols)
	    INTEGER(4), INTENT(IN) :: N, nCols
	    COMPLEX(4), INTENT(IN) :: A_mat(N,N), B_mat(N,nCols)
    END FUNCTION cSymSolve

    FUNCTION zSymSolve(A_mat, B_mat, N, nCols)
	    INTEGER(4), INTENT(IN) :: N, nCols
	    COMPLEX(8), INTENT(IN) :: A_mat(N,N), B_mat(N,nCols)
    END FUNCTION zSymSolve

    FUNCTION cHermSolve(A_mat, B_mat, N, nCols)
	    INTEGER(4), INTENT(IN) :: N, nCols
	    COMPLEX(4), INTENT(IN) :: A_mat(N,N), B_mat(N,nCols)
    END FUNCTION cHermSolve

    FUNCTION zHermSolve(A_mat, B_mat, N, nCols)
	    INTEGER(4), INTENT(IN) :: N, nCols
	    COMPLEX(8), INTENT(IN) :: A_mat(N,N), B_mat(N,nCols)
    END FUNCTION zHermSolve

    SUBROUTINE sEig(A_mat, EigVal_vec, EigVec_mat, N)
	    INTEGER(4), INTENT(IN) :: N
	    REAL(4), INTENT(IN) :: A_mat(N,N)
	    COMPLEX(4), INTENT(OUT) :: EigVal_vec(N), EigVec_mat(N,N)
    END SUBROUTINE sEig

    SUBROUTINE dEig(A_mat, EigVal_vec, EigVec_mat, N)
	    INTEGER(4), INTENT(IN) :: N
	    REAL(8), INTENT(IN) :: A_mat(N,N)
	    COMPLEX(8), INTENT(OUT) :: EigVal_vec(N), EigVec_mat(N,N)
    END SUBROUTINE dEig

    SUBROUTINE cEig(A_mat, EigVal_vec, EigVec_r_mat, N)
	    INTEGER(4), INTENT(IN) :: N
	    COMPLEX(4), INTENT(IN) :: A_mat(N,N)
	    COMPLEX(4), INTENT(OUT) :: EigVal_vec(N), EigVec_r_mat(N,N)
    END SUBROUTINE cEig

    SUBROUTINE zEig(A_mat, EigVal_vec, EigVec_r_mat, N)
	    INTEGER(4), INTENT(IN) :: N
	    COMPLEX(8), INTENT(IN) :: A_mat(N,N)
	    COMPLEX(8), INTENT(OUT) :: EigVal_vec(N), EigVec_r_mat(N,N)
    END SUBROUTINE zEig

    SUBROUTINE sSymEig(A_mat, EigVal_vec, EigVec_mat, N)
	    INTEGER(4), INTENT(IN) :: N
	    REAL(4), INTENT(IN) :: A_mat(N,N)
	    REAL(4), INTENT(OUT) :: EigVal_vec(N), EigVec_mat(N,N)
    END SUBROUTINE sSymEig

    SUBROUTINE dSymEig(A_mat, EigVal_vec, EigVec_mat, N)
	    INTEGER(4), INTENT(IN) :: N
	    REAL(8), INTENT(IN) :: A_mat(N,N)
	    REAL(8), INTENT(OUT) :: EigVal_vec(N), EigVec_mat(N,N)
    END SUBROUTINE dSymEig

    SUBROUTINE cHermEig(A_mat, EigVal_vec, EigVec_mat, N)
	    INTEGER(4), INTENT(IN) :: N
	    COMPLEX(4), INTENT(IN) :: A_mat(N,N)
	    REAL(4), INTENT(OUT) :: EigVal_vec(N)
	    COMPLEX(4), INTENT(OUT) :: EigVec_mat(N,N)
    END SUBROUTINE cHermEig

    SUBROUTINE dHermEig(A_mat, EigVal_vec, EigVec_mat, N)
	    INTEGER(4), INTENT(IN) :: N
	    COMPLEX(8), INTENT(IN) :: A_mat(N,N)
	    REAL(8), INTENT(OUT) :: EigVal_vec(N)
	    COMPLEX(8), INTENT(OUT) :: EigVec_mat(N,N)
    END SUBROUTINE dHermEig

    SUBROUTINE sGenEig(A_mat, B_mat, EigVal_vec, EigVec_mat, N)
	    INTEGER(4), INTENT(IN) :: N
	    REAL(4), INTENT(IN) :: A_mat(N,N), B_mat(N,N)
	    COMPLEX(4), INTENT(OUT) :: EigVal_vec(N), EigVec_mat(N,N)
    END SUBROUTINE sGenEig

    SUBROUTINE dGenEig(A_mat, B_mat, EigVal_vec, EigVec_mat, N)
	    INTEGER(4), INTENT(IN) :: N
	    REAL(8), INTENT(IN) :: A_mat(N,N), B_mat(N,N)
	    COMPLEX(8), INTENT(OUT) :: EigVal_vec(N), EigVec_mat(N,N)
    END SUBROUTINE dGenEig

    SUBROUTINE cGenEig(A_mat, B_mat, EigVal_vec, EigVec_r_mat, N)
	    INTEGER(4), INTENT(IN) :: N
	    COMPLEX(4), INTENT(IN) :: A_mat(N,N), B_mat(N,N)
	    COMPLEX(4), INTENT(OUT) :: EigVal_vec(N), EigVec_r_mat(N,N)
    END SUBROUTINE cGenEig

    SUBROUTINE zGenEig(A_mat, B_mat, EigVal_vec, EigVec_r_mat, N)
	    INTEGER(4), INTENT(IN) :: N
	    COMPLEX(8), INTENT(IN) :: A_mat(N,N), B_mat(N,N)
	    COMPLEX(8), INTENT(OUT) :: EigVal_vec(N), EigVec_r_mat(N,N)
    END SUBROUTINE zGenEig

    SUBROUTINE sGenSymDefEig(A_mat, B_mat, EigVal_vec, EigVec_mat, N)
	    INTEGER(4), INTENT(IN) :: N
	    REAL(4), INTENT(IN) :: A_mat(N,N), B_mat(N,N)
	    REAL(4), INTENT(OUT) :: EigVal_vec(N), EigVec_mat(N,N)
    END SUBROUTINE sGenSymDefEig

    SUBROUTINE dGenSymDefEig(A_mat, B_mat, EigVal_vec, EigVec_mat, N)
	    INTEGER(4), INTENT(IN) :: N
	    REAL(8), INTENT(IN) :: A_mat(N,N), B_mat(N,N)
	    REAL(8), INTENT(OUT) :: EigVal_vec(N), EigVec_mat(N,N)
    END SUBROUTINE dGenSymDefEig

    SUBROUTINE cGenHermDefEig(A_mat, B_mat, EigVal_vec, EigVec_mat, N)
	    INTEGER(4), INTENT(IN) :: N
	    COMPLEX(4), INTENT(IN) :: A_mat(N,N), B_mat(N,N)
	    REAL(4), INTENT(OUT) :: EigVal_vec(N)
	    COMPLEX(4), INTENT(OUT) :: EigVec_mat(N,N)
    END SUBROUTINE cGenHermDefEig

    SUBROUTINE zGenHermDefEig(A_mat, B_mat, EigVal_vec, EigVec_mat, N)
	    INTEGER(4), INTENT(IN) :: N
	    COMPLEX(8), INTENT(IN) :: A_mat(N,N), B_mat(N,N)
	    REAL(8), INTENT(OUT) :: EigVal_vec(N)
	    COMPLEX(8), INTENT(OUT) :: EigVec_mat(N,N)
    END SUBROUTINE zGenHermDefEig

END INTERFACE

END MODULE Matrix_LAPACK_interfaces