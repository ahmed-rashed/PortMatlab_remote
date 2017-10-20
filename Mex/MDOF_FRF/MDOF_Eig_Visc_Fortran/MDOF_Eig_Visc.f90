SUBROUTINE MDOF_Eig_Visc(M_mat,C_mat,K_mat,N,isProportional,EigValues_vec,EigVectors_Normalized)
    USE IFPORT
    USE Matrix_LAPACK_interfaces
    USE Matrix_Operations_interfaces

    IMPLICIT NONE
    INTEGER(4), INTENT(IN) :: N
    LOGICAL, INTENT(IN) :: isProportional
    REAL(8), INTENT(IN) ::  M_mat(N,N),C_mat(N,N),K_mat(N,N)
    COMPLEX(8), INTENT(OUT) :: EigVectors_Normalized(N,2*N), EigValues_vec(2*N)

    COMPLEX(8) i
    REAL(8), ALLOCATABLE :: EigVectors_U_mat(:,:), EigValues_U_vec(:), M_r_vec(:), w_U_r_vec(:), C_r_vec(:), zeta_r_vec(:), w_d_r_vec(:), Temp_Vec(:)
    COMPLEX(8), ALLOCATABLE :: EigValues_vec_temp(:)
    INTEGER(4) iindex(N), iiindex(2*N), NN
    INTEGER(SIZEOF_SIZE_T) array_len, array_size

    PARAMETER (i=(0D0,1D0))

    IF (isProportional .OR. all(C_mat==0)) THEN    !Undamped or proportional
        ALLOCATE(EigVectors_U_mat(N,N), EigValues_U_vec(N), M_r_vec(N), w_U_r_vec(N), C_r_vec(N), zeta_r_vec(N), w_d_r_vec(N), EigValues_vec_temp(N), Temp_Vec(N))
        CALL dGenSymDefEig(-K_mat,M_mat,EigValues_U_vec,EigVectors_U_mat,N)
    
        M_r_vec=dMatDiag(MATMUL(MATMUL(TRANSPOSE(EigVectors_U_mat),M_mat),EigVectors_U_mat),N)
        w_U_r_vec=DSQRT(-EigValues_U_vec)

        C_r_vec=dMatDiag(MATMUL(MATMUL(TRANSPOSE(EigVectors_U_mat),C_mat),EigVectors_U_mat),N)
        zeta_r_vec=C_r_vec/(2*M_r_vec*w_U_r_vec)
        w_d_r_vec=w_U_r_vec*DSQRT(1-zeta_r_vec**2)
        EigValues_vec_temp=-zeta_r_vec*w_U_r_vec+i*w_d_r_vec

        EigValues_vec(1:2*N-1:2)=EigValues_vec_temp
        EigValues_vec(2:2*N:2)=DCONJG(EigValues_vec_temp)
    
        EigVectors_Normalized(:,1:2*N-1:2)=MATMUL(EigVectors_U_mat,zDiagMat(1/zsqrt(i*2*w_d_r_vec*M_r_vec),N))
        EigVectors_Normalized(:,2:2*N:2)=DCONJG(EigVectors_Normalized(:,1:2*N-1:2))
    ELSE    !Non-proportional
        CALL quad_eig(K_mat,C_mat,M_mat,N,EigVectors_Normalized,EigValues_vec)
    END IF

    Do NN=1,N
        iindex(NN)=2*NN-1
    END DO
    array_len=N
    array_size=4
    CALL QSORT(iindex,array_len,array_size,cmp_function)
    iiindex(1:2*N-1:2)=iindex
    iiindex(2:2*N:2)=iindex+1

    EigValues_vec=EigValues_vec(iiindex)
    EigVectors_Normalized=EigVectors_Normalized(:,iiindex)

    CONTAINS
    INTEGER(2) FUNCTION cmp_function(ind1, ind2)
        INTEGER(4), INTENT(IN) :: ind1, ind2

        IF (DABS(DIMAG(EigValues_vec(ind1))) > DABS(DIMAG(EigValues_vec(ind2)))) THEN
            cmp_function=1
        ELSE IF (DABS(DIMAG(EigValues_vec(ind1))) < DABS(DIMAG(EigValues_vec(ind2)))) THEN
            cmp_function=-1
        ELSE
            cmp_function=0
        END IF

    END FUNCTION cmp_function
END SUBROUTINE MDOF_Eig_Visc