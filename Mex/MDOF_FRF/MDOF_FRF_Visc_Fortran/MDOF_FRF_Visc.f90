FUNCTION MDOF_FRF_Visc(EigValues_vec, EigVectors_Normalized_mat, N, w_column, N_w, n_vec, m_vec, N_cols)
USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_SIZE_T

IMPLICIT NONE
INTEGER(C_SIZE_T), INTENT(IN) :: N, N_w, N_cols
INTEGER(C_SIZE_T), INTENT(IN) :: n_vec(N_cols), m_vec(N_cols)
COMPLEX(8), INTENT(IN) :: EigValues_vec(2*N), EigVectors_Normalized_mat(N,2*N)
REAL(8), INTENT(IN) :: w_column(N_w,1)
COMPLEX(8) MDOF_FRF_Visc(N_w,N_cols)

INTEGER(C_SIZE_T) ii,r
COMPLEX(8) i,H_w_n_m_cols_SDOF(N_w,N_cols), A_r(2*N,2*N), A_r_temp_row(1,N_cols)

PARAMETER (i=(0D0,1D0))

    H_w_n_m_cols_SDOF=0
    MDOF_FRF_Visc=0
    DO r=1,2*N
        A_r=MATMUL(RESHAPE(EigVectors_Normalized_mat(:,r),(/N,1/)),RESHAPE(EigVectors_Normalized_mat(:,r),(/1,N/)))
        DO ii=1,N_cols
            A_r_temp_row(1,ii)=A_r(n_vec(ii), m_vec(ii))
        END DO
        H_w_n_m_cols_SDOF=H_w_n_m_cols_SDOF+MATMUL((1/(i*w_column-EigValues_vec(r))),A_r_temp_row)
    
        IF ((DIMAG(EigValues_vec(r))/=0) .AND. (MOD(r,2)/=0)) THEN !complex eigenvalue and odd r
            CYCLE
        END IF
    
        MDOF_FRF_Visc=MDOF_FRF_Visc+H_w_n_m_cols_SDOF
        H_w_n_m_cols_SDOF=0
    END DO

END FUNCTION MDOF_FRF_Visc