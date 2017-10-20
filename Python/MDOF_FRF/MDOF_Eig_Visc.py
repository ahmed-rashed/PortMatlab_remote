def MDOF_Eig_Visc(M_mat,C_mat,K_mat,isPropotional=False):
    import numpy as np
    import scipy as sp
    from quad_eig import quad_eig

    N=M_mat.shape[0]    
    if isPropotional or np.all(C_mat==0): #undamped or proportional
        EigValues_U_vec,EigVectors_U=sp.linalg.eigh(-K_mat,M_mat)    #Ordinary or generalized eigenvalue problem of a complex Hermitian or real symmetric matrix
        M_r_vec=np.diag(np.dot(np.dot(EigVectors_U.T,M_mat),EigVectors_U))
        w_U_r_vec=np.sqrt(-EigValues_U_vec)
        C_r_vec=np.diag(np.dot(np.dot(EigVectors_U.T,C_mat),EigVectors_U))
        zeta_r_vec=C_r_vec/(2*M_r_vec*w_U_r_vec)
        w_d_r_vec=w_U_r_vec*np.sqrt(1-zeta_r_vec**2)
        EigValues_vec_temp=-zeta_r_vec*w_U_r_vec+1j*w_d_r_vec
        
        EigValues_vec=np.zeros((2*N),dtype=complex)
        EigValues_vec[0:2*N-1:2]=EigValues_vec_temp
        EigValues_vec[1:2*N:2]=EigValues_vec_temp.conj()
        EigVector_Normalized=np.zeros((N,2*N),dtype=complex)
        EigVector_Normalized[:,0:2*N-1:2]=sp.linalg.solve(np.diag(np.sqrt(1j*2*w_d_r_vec*M_r_vec)),EigVectors_U.T).T
        EigVector_Normalized[:,1:2*N:2]=EigVector_Normalized[:,0:2*N-1:2].conj()
    else:    #non proportional
        EigVector_Normalized,EigValues_vec=quad_eig(K_mat,C_mat,M_mat)

    Index=np.argsort(np.absolute(EigValues_vec[0:2*N-1:2].imag))
    IIndex=np.empty(2*N,dtype=int)
    IIndex[0:2*N-1:2]=2*Index;
    IIndex[1:2*N:2]=2*Index+1;
    EigValues_vec1=EigValues_vec[IIndex]
    EigVector_Normalized1=EigVector_Normalized[:,IIndex]

    return EigVector_Normalized1,EigValues_vec1