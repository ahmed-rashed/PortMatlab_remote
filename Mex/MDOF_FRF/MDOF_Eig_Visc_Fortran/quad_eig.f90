SUBROUTINE quad_eig(K_mat,C_mat,M_mat,N,Epsi_normalized_mat,EigVal_vec)
    USE Matrix_LAPACK_interfaces
    USE Matrix_Operations_interfaces

    IMPLICIT NONE
    INTEGER(4), INTENT(IN) :: N
    REAL(8), INTENT(IN) ::  M_mat(N,N),C_mat(N,N),K_mat(N,N)
    COMPLEX(8), INTENT(OUT) :: Epsi_normalized_mat(N,2*N), EigVal_vec(2*N)

    INTEGER(4) NN
    REAL(8) N_mat(N,N),A_mat(2*N,2*N),B_mat(2*N,2*N)
    COMPLEX(8) Phi_mat(2*N,2*N), Phi_normalized_mat(2*N,2*N), B_r_vec(2*N)

    N_mat=M_mat
    A_mat=0
    A_mat(1:N,1:N)=-K_mat
    A_mat(N+1:2*N,N+1:2*N)=N_mat
    B_mat=0
    B_mat(1:N,1:N)=C_mat
    B_mat(1:N,N+1:2*N)=M_mat
    B_mat(N+1:2*N,1:N)=N_mat
    CALL dGenEig(A_mat, B_mat, EigVal_vec, Phi_mat, 2*N)
    B_r_vec=zMatDiag(MATMUL(MATMUL(TRANSPOSE(Phi_mat),B_mat),Phi_mat),2*N)

    Phi_normalized_mat=Phi_mat
    DO NN=1,2*N
        Phi_normalized_mat(:,NN)=Phi_normalized_mat(:,NN)/zsqrt(B_r_vec(NN))
    END DO

    Epsi_normalized_mat=Phi_normalized_mat(1:N,:)

END SUBROUTINE quad_eig