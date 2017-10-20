def N_DOF_sys(m_vec,c_vec,k_vec):
    import numpy as np
    
    N=m_vec.size
    
    M=np.diag(m_vec)
    C=np.diag(c_vec[0:N]+c_vec[1:N+1])+np.diag(-c_vec[1:N],1)+np.diag(-c_vec[1:N],-1)
    K=np.diag(k_vec[0:N]+k_vec[1:N+1])+np.diag(-k_vec[1:N],1)+np.diag(-k_vec[1:N],-1)
    
    return M,C,K