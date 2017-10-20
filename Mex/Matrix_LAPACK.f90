INCLUDE 'lapack.f90'

!! The following 4 FUNCTIONS call Fortran 77 LAPACK SUBROUTINES
FUNCTION sInv(A,N)
	IMPLICIT NONE
	!DEC$ NOFREEFORM
		include 'mkl_lapack.fi'
    !DEC$ FREEFORM

	INTEGER(4), INTENT(IN) :: N
	REAL(4), INTENT(IN) :: A(N,N)
	REAL(4) sInv(N,N)

	INTEGER(4) IPIV(N), LWORK, INFO
	REAL(4) L_WORK_temp(1)
	REAL(4), ALLOCATABLE :: WORK(:)

	sInv=A

	LWORK=-1
	CALL sgetri(N, sInv, N, IPIV, L_WORK_temp, LWORK, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('Lapack sgetri failed.')
	END IF

	LWORK=L_WORK_temp(1)
	ALLOCATE(WORK(LWORK))

	CALL sgetrf(N, N,sInv,N,IPIV,INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('Lapack sgetrf failed.')
	END IF

	CALL sgetri(N, sInv, N, IPIV, WORK, LWORK, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('Lapack sgetri failed.')
	END IF
END FUNCTION sInv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION dInv(A,N)
	IMPLICIT NONE
	!DEC$ NOFREEFORM
		include 'mkl_lapack.fi'
    !DEC$ FREEFORM

	INTEGER(4), INTENT(IN) :: N
	REAL(8), INTENT(IN) :: A(N,N)
	REAL(8) dInv(N,N)

	INTEGER(4) IPIV(N), LWORK, INFO
	REAL(8) L_WORK_temp(1)
	REAL(8), ALLOCATABLE :: WORK(:)

	dInv=A

	LWORK=-1
	CALL dgetri(N, dInv, N, IPIV, L_WORK_temp, LWORK, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('Lapack dgetri failed.')
	END IF

	LWORK=L_WORK_temp(1)
	ALLOCATE(WORK(LWORK))

	CALL dgetrf(N, N,dInv,N,IPIV,INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('Lapack dgetrf failed.')
	END IF

	CALL dgetri(N, dInv, N, IPIV, WORK, LWORK, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('Lapack dgetri failed.')
	END IF
END FUNCTION dInv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION cInv(A,N)
	IMPLICIT NONE
	!DEC$ NOFREEFORM
		include 'mkl_lapack.fi'
    !DEC$ FREEFORM

	INTEGER(4), INTENT(IN) :: N
	COMPLEX(4), INTENT(IN) :: A(N,N)
	COMPLEX(4) cInv(N,N)

	INTEGER(4) IPIV(N), LWORK, INFO
	COMPLEX(4) L_WORK_temp(1)
	COMPLEX(4), ALLOCATABLE :: WORK(:)

	cInv=A

	LWORK=-1
	CALL cgetri(N, cInv, N, IPIV, L_WORK_temp, LWORK, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('Lapack cgetri failed.')
	END IF

	LWORK=REAL(L_WORK_temp(1))
	ALLOCATE(WORK(LWORK))

	CALL cgetrf(N, N,cInv,N,IPIV,INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('Lapack cgetrf failed.')
	END IF

	CALL cgetri(N, cInv, N, IPIV, WORK, LWORK, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('Lapack cgetri failed.')
	END IF
END FUNCTION cInv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION zInv(A,N)
	IMPLICIT NONE
	!DEC$ NOFREEFORM
		include 'mkl_lapack.fi'
    !DEC$ FREEFORM

	INTEGER(4), INTENT(IN) :: N
	COMPLEX(8), INTENT(IN) :: A(N,N)
	COMPLEX(8) zInv(N,N)

	INTEGER(4) IPIV(N), LWORK, INFO
	REAL(8) L_WORK_temp(1)
	REAL(8), ALLOCATABLE :: WORK(:)

	zInv=A

	LWORK=-1
	CALL zgetri(N, zInv, N, IPIV, L_WORK_temp, LWORK, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('Lapack zgetri failed.')
	END IF

	LWORK=L_WORK_temp(1)
	ALLOCATE(WORK(LWORK))

	CALL zgetrf(N, N,zInv,N,IPIV,INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('Lapack zgetrf failed.')
	END IF

	CALL zgetri(N, zInv, N, IPIV, WORK, LWORK, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('Lapack zgetri failed.')
	END IF
END FUNCTION zInv

!! All the following function/Subroutines use Fortran 95 LAPACK subroutines
FUNCTION sSymInv(A_mat, N)
	USE lapack95, ONLY: sytrf, sytri
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N
	REAL(4), INTENT(IN) :: A_mat(N,N)
	REAL(4) sSymInv(N,N)
	INTEGER(4) ipiv(N), INFO

	sSymInv=A_mat
	CALL sytrf(sSymInv, 'U' , ipiv, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('sytrf failed.')
	END IF

	CALL sytri(sSymInv, ipiv, 'U', INFO)
	IF (INFO /= 0) THEN
	    CALL mexErrMsgTxt('sytri failed.')
	END IF
END FUNCTION sSymInv

FUNCTION dSymInv(A_mat, N)
	USE lapack95, ONLY: sytrf, sytri
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N
	REAL(8), INTENT(IN) :: A_mat(N,N)
	REAL(8) dSymInv(N,N)
	INTEGER(4) ipiv(N), INFO
	
	dSymInv=A_mat
	CALL sytrf(dSymInv, 'U', ipiv, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('sytrf failed.')
	END IF
	
	CALL sytri(dSymInv, ipiv, 'U', INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('sytri failed.')
	END IF
END FUNCTION dSymInv

FUNCTION cSymInv(A_mat, N)
	USE lapack95, ONLY: sytrf, sytri
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N
	COMPLEX(4), INTENT(IN) :: A_mat(N,N)
	COMPLEX(4) cSymInv(N,N)
	INTEGER(4) ipiv(N), INFO

	cSymInv=A_mat
	CALL sytrf(cSymInv, 'U', ipiv, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('sytrf failed.')
	END IF
	
	CALL sytri(cSymInv, ipiv, 'U', INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('sytri failed.')
	END IF

END FUNCTION cSymInv

FUNCTION zSymInv(A_mat, N)
	USE lapack95, ONLY: sytrf, sytri
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N
	COMPLEX(8), INTENT(IN) :: A_mat(N,N)
	COMPLEX(8) zSymInv(N,N)
	INTEGER(4) ipiv(N), INFO
	
	zSymInv=A_mat
	CALL sytrf(zSymInv, 'U', ipiv, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('sytrf failed.')
	END IF
	
	CALL sytri(zSymInv, ipiv, 'U', INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('sytri failed.')
	END IF
END FUNCTION zSymInv

FUNCTION sSymPosDefInv(A_mat, N)
	USE lapack95, ONLY: potrf, potri
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N
	REAL(4), INTENT(IN) :: A_mat(N,N)
	REAL(4) sSymPosDefInv(N,N)
    INTEGER(4) INFO

	sSymPosDefInv=A_mat
	CALL potrf(sSymPosDefInv, 'U', INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('potrf failed.')
	END IF
	CALL potri(sSymPosDefInv, 'U', INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('potri failed.')
	END IF

END FUNCTION sSymPosDefInv

FUNCTION dSymPosDefInv(A_mat, N)
	USE lapack95, ONLY: potrf, potri
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N
	REAL(8), INTENT(IN) :: A_mat(N,N)
	REAL(8) dSymPosDefInv(N,N)
    INTEGER(4) INFO
	dSymPosDefInv=A_mat
	CALL potrf(dSymPosDefInv, 'U', INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('potrf failed.')
	END IF
	
	CALL potri(dSymPosDefInv, 'U', INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('potri failed.')
	END IF
END FUNCTION dSymPosDefInv

FUNCTION cHermPosDefInv(A_mat, N)
	USE lapack95, ONLY: potrf, potri
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N
	COMPLEX(4), INTENT(IN) :: A_mat(N,N)
	COMPLEX(4) cHermPosDefInv(N,N)
    INTEGER(4) INFO

	cHermPosDefInv=A_mat
	CALL potrf(cHermPosDefInv, 'U', INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('potrf failed.')
	END IF
	
	CALL potri(cHermPosDefInv, 'U', INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('potri failed.')
	END IF
END FUNCTION cHermPosDefInv

FUNCTION zHermPosDefInv(A_mat, N)
	USE lapack95, ONLY: potrf, potri
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N
	COMPLEX(8), INTENT(IN) :: A_mat(N,N)
	COMPLEX(8) zHermPosDefInv(N,N)
    INTEGER(4) INFO

	zHermPosDefInv=A_mat
	CALL potrf(zHermPosDefInv, 'U', INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('potrf failed.') 
	END IF
	
	CALL potri(zHermPosDefInv, 'U', INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('potri failed.') 
	END IF
END FUNCTION zHermPosDefInv


FUNCTION sSolve(A_mat, B_mat, N, nCols)
	USE lapack95, ONLY: gesv
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N, nCols
	REAL(4), INTENT(IN) :: A_mat(N,N), B_mat(N,nCols)
	REAL(4) sSolve(N,nCols)
	REAL(4) A_temp(N,N)
	INTEGER(4) ipiv(N), INFO

	A_temp=A_mat
	sSolve=B_mat
	CALL gesv(A_temp, sSolve, ipiv, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('gesv failed.') 
	END IF
END FUNCTION sSolve

FUNCTION dSolve(A_mat, B_mat, N, nCols)
	USE lapack95, ONLY: gesv
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N, nCols
	REAL(8), INTENT(IN) :: A_mat(N,N), B_mat(N,nCols)
	REAL(8) dSolve(N,nCols)
	REAL(8) A_temp(N,N)
	INTEGER(4) ipiv(N), INFO

	A_temp=A_mat
	dSolve=B_mat
	CALL gesv(A_temp, dSolve, ipiv, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('gesv failed.') 
	END IF
END FUNCTION dSolve

FUNCTION cSolve(A_mat, B_mat, N, nCols)
	USE lapack95, ONLY: gesv
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N, nCols
	COMPLEX(4), INTENT(IN) :: A_mat(N,N), B_mat(N,nCols)
	COMPLEX(4) cSolve(N,nCols)
	COMPLEX(4) A_temp(N,N)
	INTEGER(4) ipiv(N), INFO
	
	A_temp=A_mat
	cSolve=B_mat
	CALL gesv(A_temp, cSolve, ipiv, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('gesv failed.') 
	END IF
END FUNCTION cSolve

FUNCTION zSolve(A_mat, B_mat, N, nCols)
	USE lapack95, ONLY: gesv
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N, nCols
	COMPLEX(8), INTENT(IN) :: A_mat(N,N), B_mat(N,nCols)
	COMPLEX(8) zSolve(N,nCols)
	COMPLEX(8) A_temp(N,N)
	INTEGER(4) ipiv(N), INFO
	
	A_temp=A_mat
	zSolve=B_mat
	CALL gesv(A_temp, zSolve, ipiv, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('gesv failed.') 
	END IF
END FUNCTION zSolve

FUNCTION sSymPosDefSolve(A_mat, B_mat, N, nCols)
	USE lapack95, ONLY: posv
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N, nCols
	REAL(4), INTENT(IN) :: A_mat(N,N), B_mat(N,nCols)
	REAL(4) sSymPosDefSolve(N,nCols)
	REAL(4) A_temp(N,N)
	INTEGER(4) INFO
	A_temp=A_mat
	sSymPosDefSolve=B_mat
	CALL posv(A_temp, sSymPosDefSolve, 'U', INFO)
	IF (INFO > 0) THEN
		CALL mexErrMsgTxt('Error using posv, matrix A_mat is not positive definite.')
	ELSE IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('posv failed.')
	END IF
END FUNCTION sSymPosDefSolve

FUNCTION dSymPosDefSolve(A_mat, B_mat, N, nCols)
	USE lapack95, ONLY: posv
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N, nCols
	REAL(8), INTENT(IN) :: A_mat(N,N), B_mat(N,nCols)
	REAL(8) dSymPosDefSolve(N,nCols)
	REAL(8) A_temp(N,N)
    INTEGER(4) INFO
	
	A_temp=A_mat
	dSymPosDefSolve=B_mat
	CALL posv(A_temp, dSymPosDefSolve, 'U', INFO)
	IF (INFO > 0) THEN
		CALL mexErrMsgTxt('Error using posv, matrix A_mat is not positive definite.')
	ELSE IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('posv failed.')
	END IF
END FUNCTION dSymPosDefSolve

FUNCTION cHermPosDefSolve(A_mat, B_mat, N, nCols)
	USE lapack95, ONLY: posv
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N, nCols
	COMPLEX(4), INTENT(IN) :: A_mat(N,N), B_mat(N,nCols)
	COMPLEX(4) cHermPosDefSolve(N,nCols)
	COMPLEX(4) A_temp(N,N)
    INTEGER(4) INFO
	
	A_temp=A_mat
	cHermPosDefSolve=B_mat
	CALL posv(A_temp, cHermPosDefSolve, 'U', INFO)
	IF (INFO > 0) THEN
		CALL mexErrMsgTxt('Error using posv, matrix A_mat is not positive definite.')
	ELSE IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('posv failed.')
	END IF
END FUNCTION cHermPosDefSolve

FUNCTION zHermPosDefSolve(A_mat, B_mat, N, nCols)
	USE lapack95, ONLY: posv
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N, nCols
	COMPLEX(8), INTENT(IN) :: A_mat(N,N), B_mat(N,nCols)
	COMPLEX(8) zHermPosDefSolve(N,nCols)
	COMPLEX(8) A_temp(N,N)
    INTEGER(4) INFO
	
	A_temp=A_mat
	zHermPosDefSolve=B_mat
	CALL posv(A_temp, zHermPosDefSolve, 'U', INFO)
	IF (INFO > 0) THEN
		CALL mexErrMsgTxt('Error using posv, matrix A_mat is not positive definite.')
	ELSE IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('posv failed.')
	END IF
END FUNCTION zHermPosDefSolve

FUNCTION sSymSolve(A_mat, B_mat, N, nCols)
	USE lapack95, ONLY: sysv
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N, nCols
	REAL(4), INTENT(IN) :: A_mat(N,N), B_mat(N,nCols)
	REAL(4) sSymSolve(N,nCols)
	INTEGER(4) ipiv(N)
	REAL(4) A_temp(N,N)
    INTEGER(4) INFO
	
	A_temp=A_mat
	sSymSolve=B_mat
	CALL sysv(A_temp, sSymSolve, 'U', ipiv, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('sysv failed.') 
	END IF
END FUNCTION sSymSolve

FUNCTION dSymSolve(A_mat, B_mat, N, nCols)
	USE lapack95, ONLY: sysv
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N, nCols
	REAL(8), INTENT(IN) :: A_mat(N,N), B_mat(N,nCols)
	REAL(8) dSymSolve(N,nCols)
	INTEGER(4) ipiv(N)
	REAL(8) A_temp(N,N)
    INTEGER(4) INFO
	
	A_temp=A_mat
	dSymSolve=B_mat
	CALL sysv(A_temp, dSymSolve, 'U', ipiv, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('sysv failed.') 
	END IF
END FUNCTION dSymSolve

FUNCTION cSymSolve(A_mat, B_mat, N, nCols)
	USE lapack95, ONLY: sysv
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N, nCols
	COMPLEX(4), INTENT(IN) :: A_mat(N,N), B_mat(N,nCols)
	COMPLEX(4) cSymSolve(N,nCols)
	INTEGER(4) ipiv(N)
	COMPLEX(4) A_temp(N,N)
	INTEGER(4) INFO

	A_temp=A_mat
	cSymSolve=B_mat
	CALL sysv(A_temp, cSymSolve, 'U', ipiv, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('sysv failed.') 
	END IF
END FUNCTION cSymSolve

FUNCTION zSymSolve(A_mat, B_mat, N, nCols)
	USE lapack95, ONLY: sysv
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N, nCols
	COMPLEX(8), INTENT(IN) :: A_mat(N,N), B_mat(N,nCols)
	COMPLEX(8) zSymSolve(N,nCols)
	INTEGER(4) ipiv(N)
	COMPLEX(8) A_temp(N,N)
    INTEGER(4) INFO
	
	A_temp=A_mat
	zSymSolve=B_mat
	CALL sysv(A_temp, zSymSolve, 'U', ipiv, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('sysv failed.') 
	END IF
END FUNCTION zSymSolve

FUNCTION cHermSolve(A_mat, B_mat, N, nCols)
	USE lapack95, ONLY: hesv
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N, nCols
	COMPLEX(4), INTENT(IN) :: A_mat(N,N), B_mat(N,nCols)
	COMPLEX(4) cHermSolve(N,nCols)
	INTEGER(4) ipiv(N), INFO
	COMPLEX(4) A_temp(N,N)
	
	A_temp=A_mat
	cHermSolve=B_mat
	CALL hesv(A_temp, cHermSolve, 'U', ipiv, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('hesv failed.') 
	END IF
END FUNCTION cHermSolve

FUNCTION zHermSolve(A_mat, B_mat, N, nCols)
	USE lapack95, ONLY: hesv
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N, nCols
	COMPLEX(8), INTENT(IN) :: A_mat(N,N), B_mat(N,nCols)
	COMPLEX(8) zHermSolve(N,nCols)
	INTEGER(4) ipiv(N), INFO
	COMPLEX(8) A_temp(N,N)
	
	A_temp=A_mat
	zHermSolve=B_mat
	CALL hesv(A_temp, zHermSolve, 'U', ipiv, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('hesv failed.') 
	END IF
END FUNCTION zHermSolve

SUBROUTINE sEig(A_mat, EigVal_vec, EigVec_mat, N)
	USE lapack95, ONLY: geev
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N
	REAL(4), INTENT(IN) :: A_mat(N,N)
	COMPLEX(4), INTENT(OUT) :: EigVal_vec(N), EigVec_mat(N,N)

	REAL(4) A_temp(N,N), EigVec_l(N,N), EigVec_r(N,N)
	REAL(4) w_r(N), w_i(N)
	INTEGER(4) NN, INFO
	COMPLEX(4) i

    PARAMETER (i=(0D0,1D0))

    A_temp=A_mat
	CALL geev(A_temp, w_r, w_i, EigVec_l, EigVec_r, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('geev failed.') 
	END IF

	EigVal_vec=w_r+ w_i*i
	DO NN=1,N
		IF (w_i(NN)/=0) THEN
            IF (MOD(NN,2)==0) THEN
                CYCLE
            END IF
			EigVec_mat(:,NN)=EigVec_r(:,NN)+i*EigVec_r(:,NN+1)
			EigVec_mat(:,NN+1)=EigVec_r(:,NN)-i*EigVec_r(:,NN+1)
		ELSE
			EigVec_mat(:,NN)=EigVec_r(:,NN)
		END IF
	END DO
END SUBROUTINE sEig

SUBROUTINE dEig(A_mat, EigVal_vec, EigVec_mat, N)
	USE lapack95, ONLY: geev
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N
	REAL(8), INTENT(IN) :: A_mat(N,N)
	COMPLEX(8), INTENT(OUT) :: EigVal_vec(N), EigVec_mat(N,N)

	REAL(8) A_temp(N,N), EigVec_l(N,N), EigVec_r(N,N)
	REAL(8) w_r(N), w_i(N)
	INTEGER(4) NN, INFO
    COMPLEX(8) i

    PARAMETER (i=(0D0,1D0))

    A_temp=A_mat
	CALL geev(A_temp, w_r, w_i, EigVec_l, EigVec_r, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('geev failed.') 
	END IF

	EigVal_vec=w_r+ w_i*i
	DO NN=1,N
		IF (w_i(NN)/=0) THEN
            IF (MOD(NN,2)==0) THEN
                CYCLE
            END IF
			EigVec_mat(:,NN)=EigVec_r(:,NN)+i*EigVec_r(:,NN+1)
			EigVec_mat(:,NN+1)=EigVec_r(:,NN)-i*EigVec_r(:,NN+1)
		ELSE
			EigVec_mat(:,NN)=EigVec_r(:,NN)
		END IF
	END DO
END SUBROUTINE dEig

SUBROUTINE cEig(A_mat, EigVal_vec, EigVec_r_mat, N)
	USE lapack95, ONLY: geev
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N
	COMPLEX(4), INTENT(IN) :: A_mat(N,N)
	COMPLEX(4), INTENT(OUT) :: EigVal_vec(N), EigVec_r_mat(N,N)

	COMPLEX(4) A_temp(N,N), EigVec_l_mat(N,N), i
    INTEGER(4) INFO

    PARAMETER (i=(0D0,1D0))

    A_temp=A_mat
	CALL geev(A_temp, EigVal_vec, EigVec_l_mat, EigVec_r_mat, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('geev failed.') 
	END IF
END SUBROUTINE cEig

SUBROUTINE zEig(A_mat, EigVal_vec, EigVec_r_mat, N)
	USE lapack95, ONLY: geev
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N
	COMPLEX(8), INTENT(IN) :: A_mat(N,N)
	COMPLEX(8), INTENT(OUT) :: EigVal_vec(N), EigVec_r_mat(N,N)

	COMPLEX(8) A_temp(N,N), EigVec_l_mat(N,N), i
    INTEGER(4) INFO

    PARAMETER (i=(0D0,1D0))

    A_temp=A_mat
	CALL geev(A_temp, EigVal_vec, EigVec_l_mat, EigVec_r_mat, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('geev failed.') 
	END IF
END SUBROUTINE zEig

SUBROUTINE sSymEig(A_mat, EigVal_vec, EigVec_mat, N)
	USE lapack95, ONLY: syevd
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N
	REAL(4), INTENT(IN) :: A_mat(N,N)
	REAL(4), INTENT(OUT) :: EigVal_vec(N), EigVec_mat(N,N)
    INTEGER(4) INFO

	EigVec_mat=A_mat
	CALL syevd(EigVec_mat, EigVal_vec, 'V', 'U', INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('syevd failed.') 
	END IF
END SUBROUTINE sSymEig

SUBROUTINE dSymEig(A_mat, EigVal_vec, EigVec_mat, N)
	USE lapack95, ONLY: syevd
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N
	REAL(8), INTENT(IN) :: A_mat(N,N)
	REAL(8), INTENT(OUT) :: EigVal_vec(N), EigVec_mat(N,N)
    INTEGER(4) INFO

	EigVec_mat=A_mat

	CALL syevd(EigVec_mat, EigVal_vec, 'V', 'U', INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('syevd failed.') 
	END IF
END SUBROUTINE dSymEig

SUBROUTINE cHermEig(A_mat, EigVal_vec, EigVec_mat, N)
	USE lapack95, ONLY: heevd
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N
	COMPLEX(4), INTENT(IN) :: A_mat(N,N)
	REAL(4), INTENT(OUT) :: EigVal_vec(N)
	COMPLEX(4), INTENT(OUT) :: EigVec_mat(N,N)
    INTEGER(4) INFO

	EigVec_mat=A_mat
	CALL heevd(EigVec_mat, EigVal_vec, 'V', 'U', INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('heevd failed.') 
	END IF
END SUBROUTINE cHermEig

SUBROUTINE dHermEig(A_mat, EigVal_vec, EigVec_mat, N)
	USE lapack95, ONLY: heevd
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N
	COMPLEX(8), INTENT(IN) :: A_mat(N,N)
	REAL(8), INTENT(OUT) :: EigVal_vec(N)
	COMPLEX(8), INTENT(OUT) :: EigVec_mat(N,N)
    INTEGER(4) INFO

	EigVec_mat=A_mat
	CALL heevd(EigVec_mat, EigVal_vec, 'V', 'U', INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('heevd failed.') 
	END IF
END SUBROUTINE dHermEig

SUBROUTINE sGenEig(A_mat, B_mat, EigVal_vec, EigVec_mat, N)
	USE lapack95, ONLY: ggev
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N
	REAL(4), INTENT(IN) :: A_mat(N,N), B_mat(N,N)
	COMPLEX(4), INTENT(OUT) :: EigVal_vec(N), EigVec_mat(N,N)

	REAL(4) A_temp(N,N), B_temp(N,N), EigVec_l(N,N), EigVec_r(N,N)
	REAL(4) alphar(N), alphai(N), beta(N)
	INTEGER(4) NN, INFO
	COMPLEX(4) i

    PARAMETER (i=(0D0,1D0))

    A_temp=A_mat
    B_temp=B_mat
	CALL ggev(A_temp, B_temp, alphar, alphai, beta, EigVec_l, EigVec_r, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('ggev failed.') 
	END IF

	EigVal_vec=(alphar+ alphai*i)/beta
	DO NN=1,N
		IF (alphai(NN)/=0) THEN
            IF (MOD(NN,2)==0) THEN
                CYCLE
            END IF
			EigVec_mat(:,NN)=EigVec_r(:,NN)+i*EigVec_r(:,NN+1)
			EigVec_mat(:,NN+1)=EigVec_r(:,NN)-i*EigVec_r(:,NN+1)
		ELSE
			EigVec_mat(:,NN)=EigVec_r(:,NN)
		END IF
	END DO
END SUBROUTINE sGenEig

SUBROUTINE dGenEig(A_mat, B_mat, EigVal_vec, EigVec_mat, N)
	USE lapack95, ONLY: ggev
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N
	REAL(8), INTENT(IN) :: A_mat(N,N), B_mat(N,N)
	COMPLEX(8), INTENT(OUT) :: EigVal_vec(N), EigVec_mat(N,N)

	REAL(8) A_temp(N,N), B_temp(N,N), EigVec_l(N,N), EigVec_r(N,N)
	REAL(8) alphar(N), alphai(N), beta(N)
	INTEGER(4) NN, INFO
	COMPLEX(8) i

    PARAMETER (i=(0D0,1D0))

    A_temp=A_mat
    B_temp=B_mat
	CALL ggev(A_temp, B_temp, alphar, alphai, beta, EigVec_l, EigVec_r, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('ggev failed.') 
	END IF

	EigVal_vec=(alphar+alphai*i)/beta
	DO NN=1,N
		IF (alphai(NN)/=0) THEN
            IF (MOD(NN,2)==0) THEN
                CYCLE
            END IF
			EigVec_mat(:,NN)=EigVec_r(:,NN)+i*EigVec_r(:,NN+1)
			EigVec_mat(:,NN+1)=EigVec_r(:,NN)-i*EigVec_r(:,NN+1)
		ELSE
			EigVec_mat(:,NN)=EigVec_r(:,NN)
		END IF
	END DO
END SUBROUTINE dGenEig

SUBROUTINE cGenEig(A_mat, B_mat, EigVal_vec, EigVec_r_mat, N)
	USE lapack95, ONLY: ggev
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N
	COMPLEX(4), INTENT(IN) :: A_mat(N,N), B_mat(N,N)
	COMPLEX(4), INTENT(OUT) :: EigVal_vec(N), EigVec_r_mat(N,N)

	COMPLEX(4) A_temp(N,N), B_temp(N,N), EigVec_l_mat(N,N), alpha(N), beta(N), i
    INTEGER(4) INFO

    PARAMETER (i=(0D0,1D0))

    A_temp=A_mat
    B_temp=B_mat
	CALL ggev(A_temp, B_temp, alpha, beta, EigVec_l_mat, EigVec_r_mat, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('ggev failed.') 
	END IF

	EigVal_vec=alpha/beta
END SUBROUTINE cGenEig

SUBROUTINE zGenEig(A_mat, B_mat, EigVal_vec, EigVec_r_mat, N)
	USE lapack95, ONLY: ggev
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N
	COMPLEX(8), INTENT(IN) :: A_mat(N,N), B_mat(N,N)
	COMPLEX(8), INTENT(OUT) :: EigVal_vec(N), EigVec_r_mat(N,N)

	COMPLEX(8) A_temp(N,N), B_temp(N,N), EigVec_l_mat(N,N), alpha(N), beta(N), i
    INTEGER(4) INFO

    PARAMETER (i=(0D0,1D0))

    A_temp=A_mat
    B_temp=B_mat
	CALL ggev(A_temp, B_temp, alpha, beta, EigVec_l_mat, EigVec_r_mat, INFO)
	IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('ggev failed.') 
	END IF

	EigVal_vec=alpha/beta
END SUBROUTINE zGenEig

SUBROUTINE sGenSymDefEig(A_mat, B_mat, EigVal_vec, EigVec_mat, N)
	USE lapack95, ONLY: sygvd
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N
	REAL(4), INTENT(IN) :: A_mat(N,N), B_mat(N,N)
	REAL(4), INTENT(OUT) :: EigVal_vec(N), EigVec_mat(N,N)

	REAL(4) B_temp(N,N)
    INTEGER(4) INFO

	EigVec_mat=A_mat
    B_temp=B_mat
	CALL sygvd(EigVec_mat, B_temp, EigVal_vec, 1, 'V', 'U', INFO)
	IF (INFO > N) THEN
		CALL mexErrMsgTxt('Error using sygvd, matrix B_mat is not positive definite.')
	ELSE IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('sygvd failed.')
	END IF
END SUBROUTINE sGenSymDefEig

SUBROUTINE dGenSymDefEig(A_mat, B_mat, EigVal_vec, EigVec_mat, N)
	USE lapack95, ONLY: sygvd
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N
	REAL(8), INTENT(IN) :: A_mat(N,N), B_mat(N,N)
	REAL(8), INTENT(OUT) :: EigVal_vec(N), EigVec_mat(N,N)

	REAL(8) B_temp(N,N)
    INTEGER(4) INFO

	EigVec_mat=A_mat
    B_temp=B_mat
	CALL sygvd(EigVec_mat, B_temp, EigVal_vec, 1, 'V', 'U', INFO)
	IF (INFO > N) THEN
		CALL mexErrMsgTxt('Error using sygvd, matrix B_mat is not positive definite.')
	ELSE IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('sygvd failed.')
	END IF
END SUBROUTINE dGenSymDefEig

SUBROUTINE cGenHermDefEig(A_mat, B_mat, EigVal_vec, EigVec_mat, N)
	USE lapack95, ONLY: hegvd
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N
	COMPLEX(4), INTENT(IN) :: A_mat(N,N), B_mat(N,N)
	REAL(4), INTENT(OUT) :: EigVal_vec(N)
	COMPLEX(4), INTENT(OUT) :: EigVec_mat(N,N)

	COMPLEX(4) B_temp(N,N)
    INTEGER(4) INFO

	EigVec_mat=A_mat
    B_temp=B_mat
	CALL hegvd(EigVec_mat, B_temp, EigVal_vec, 1, 'V', 'U', INFO)
	IF (INFO > N) THEN
		CALL mexErrMsgTxt('Error using hegvd, matrix B_mat is not positive definite.')
	ELSE IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('hegvd failed.')
	END IF
END SUBROUTINE cGenHermDefEig

SUBROUTINE zGenHermDefEig(A_mat, B_mat, EigVal_vec, EigVec_mat, N)
	USE lapack95, ONLY: hegvd
	IMPLICIT NONE

	INTEGER(4), INTENT(IN) :: N
	COMPLEX(8), INTENT(IN) :: A_mat(N,N), B_mat(N,N)
	REAL(8), INTENT(OUT) :: EigVal_vec(N)
	COMPLEX(8), INTENT(OUT) :: EigVec_mat(N,N)

	COMPLEX(8) B_temp(N,N)
    INTEGER(4) INFO

	EigVec_mat=A_mat
    B_temp=B_mat
	CALL hegvd(EigVec_mat, B_temp, EigVal_vec, 1, 'V', 'U', INFO)
	IF (INFO > N) THEN
		CALL mexErrMsgTxt('Error using hegvd, matrix B_mat is not positive definite.')
	ELSE IF (INFO /= 0) THEN
		CALL mexErrMsgTxt('hegvd failed.')
	END IF
END SUBROUTINE zGenHermDefEig