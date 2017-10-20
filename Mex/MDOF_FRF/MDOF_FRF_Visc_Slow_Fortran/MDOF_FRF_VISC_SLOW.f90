FUNCTION MDOF_FRF_VISC_SLOW(M_mat,C_mat,K_mat,N,w_col,n_w_points,n_vec,m_vec,n_n_row)

USE Matrix_LAPACK_interfaces

IMPLICIT NONE
INTEGER(4), INTENT(IN) :: N,n_w_points,n_n_row
REAL(8), INTENT(IN) ::  M_mat(N,N),C_mat(N,N),K_mat(N,N),w_col(n_w_points,1),n_vec(n_n_row),m_vec(n_n_row)
COMPLEX(8) MDOF_FRF_VISC_SLOW(n_w_points,n_n_row)

INTEGER(4) nn,mm
COMPLEX(8) i,H_w_mat(N,N)
REAL(8) pi

PARAMETER (pi=DATAN(1D0)*4D0, i=(0D0,1D0))

	DO nn=1,n_w_points
		H_w_mat=zInv(-w_col(nn,1)*w_col(nn,1)*M_mat+i*w_col(nn,1)*C_mat+K_mat,N)

        DO mm=1,n_n_row
		    MDOF_FRF_VISC_SLOW(nn,mm)=H_w_mat(n_vec(mm),m_vec(mm))
        END DO
	END DO

END FUNCTION MDOF_FRF_VISC_SLOW