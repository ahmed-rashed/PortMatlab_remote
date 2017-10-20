def quad_eig(K_mat,C_mat,M_mat):
    import numpy as np
    import scipy as sp
    import scipy.linalg as la
    #from scipy import linalg as la
    
    n=M_mat.shape[0]
    N_mat=M_mat
    A=np.concatenate((np.concatenate((-K_mat,0*M_mat),axis=1),np.concatenate((0*M_mat,N_mat),axis=1)),axis=0)
    B=np.concatenate((np.concatenate((C_mat,M_mat),axis=1),np.concatenate((N_mat,0*M_mat),axis=1)),axis=0)
    Val_vec,Phi=la.eig(A,B) #Ordinary or generalized eigenvalue problem of a complex Hermitian or real symmetric matrix
    B_r=np.dot(np.dot(Phi.T,B),Phi)
    Phi_normalized=la.lstsq(np.sqrt(B_r.T),Phi.T)[0].T
    Epsi_normalized=Phi_normalized[0:n,:]
    
    return Epsi_normalized,Val_vec