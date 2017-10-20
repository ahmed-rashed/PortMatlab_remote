def MDOF_FRF_Visc_slow(M_mat,C_mat,K_mat,w_vec,n_row,m_row):
    import numpy as np

    n_w_points=w_vec.size
    H_cols=np.zeros((n_w_points,m_row.size),dtype=complex)
    for ii in np.arange(0,n_w_points):
        H_w=np.linalg.inv(-w_vec[ii]**2*M_mat+1j*w_vec[ii]*C_mat+K_mat)

        H_cols[ii,:]=H_w[n_row,m_row]

    return H_cols