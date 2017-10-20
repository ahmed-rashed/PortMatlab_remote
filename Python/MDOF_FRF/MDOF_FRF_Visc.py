"Initialized by Alyaa Mansour, 1st year Aerospace; Spring 2016"

def MDOF_FRF_Visc(EigValues_vec,EigVectors_Normalized,w_col,ii_vec,jj_vec):
    import numpy as np
	
    N=EigVectors_Normalized.shape[0]
    n_w=w_col.size
    if ii_vec.size!=jj_vec.size:
        raise NameError('Dimensions of ii_vec and jj_vec must be identical')

    n_col=ii_vec.size
    H_w_mat_SDOF=np.zeros((n_w,n_col),dtype=complex)
    H_w_mat=np.zeros((n_w,n_col),dtype=complex)
    for r in np.arange(2*N):
        A_r=np.dot(EigVectors_Normalized[:,[r]],EigVectors_Normalized[:,[r]].T)
        A_r_temp_row=np.reshape(A_r[ii_vec,jj_vec],(1,n_col))
        H_w_mat_SDOF=H_w_mat_SDOF+np.dot((1/(1j*w_col-EigValues_vec[r])),A_r_temp_row)
        if EigValues_vec[r].imag!=0 and np.mod(r+1,2)!=0:   #complex eigenvalue and odd r
            continue

        H_w_mat=H_w_mat+H_w_mat_SDOF
        H_w_mat_SDOF=np.zeros((n_w,n_col),dtype=complex)
        
    return H_w_mat