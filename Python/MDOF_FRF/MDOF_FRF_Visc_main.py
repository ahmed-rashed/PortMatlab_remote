import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from MDOF_FRF_Visc_slow import MDOF_FRF_Visc_slow
from MDOF_FRF_Visc import MDOF_FRF_Visc
from MDOF_Eig_Visc import MDOF_Eig_Visc
from N_DOF_sys import N_DOF_sys
from mpl_toolkits.mplot3d import Axes3D
from plot_FRF_3d import plot_FRF_3d
from plot_FRF_Nyq import plot_FRF_Nyq
from plot_FRF_r_i import plot_FRF_r_i
from plot_FRF_mag_phase import plot_FRF_mag_phase

if __name__ == '__main__':
    plt.close('all')
    matplotlib.rcParams['text.usetex'] = True
    
    f_final=40
    n_f_points=8000
    
    N=3
    m_vec=10*np.ones(N)
    c_vec=8*np.ones(N+1)
    k_vec=np.concatenate((161000/2*np.ones(N),[0]), axis=0)
    
    n_row=np.array([0,0,0])
    m_row=np.array([2,1,0])
    
    isproportional=False
    M_mat,C_mat,K_mat=N_DOF_sys(m_vec,c_vec,k_vec)
    #C_mat=0;isproportional=True
    #C_mat=.0001*M_mat+.0001*K_mat;isproportional=True
    
    f_col=np.linspace(0,f_final,n_f_points).reshape(n_f_points,1)
    w_col=2*np.pi*f_col

    #Slow FRF calculation
    H_cols0=MDOF_FRF_Visc_slow(M_mat,C_mat,K_mat,w_col,n_row,m_row)

    #Fast FRF calculation
    EigVectors_Normalized, EigValues_vec=MDOF_Eig_Visc(M_mat, C_mat, K_mat,isproportional)
    H_cols=MDOF_FRF_Visc(EigValues_vec, EigVectors_Normalized, w_col, n_row, m_row)

    print np.amax(np.absolute(H_cols-H_cols0))

    f_3d=plt.figure()
    f_3d.gca(projection='3d')
    
    f_Nyq=plt.figure()
    f_r_i=plt.figure()
    f_mag_phase=plt.figure()
    ax_mag=None
    ax_phase=None
    for ii in range(m_row.size):
        plt.figure(f_3d.number)
        plot_FRF_3d(f_col,H_cols[:,ii],lineLabel=r'$H_{%u,%u}$' %(n_row[ii]+1,m_row[ii]+1),DispRealImag=True)[0]

        plt.figure(f_Nyq.number)
        plot_FRF_Nyq(H_cols[:,ii],lineLabel=r'$H_{%u,%u}$' %(n_row[ii]+1,m_row[ii]+1))[0]
        
        plt.figure(f_r_i.number)
        ax_r,ax_i=plot_FRF_r_i(f_col,H_cols[:,ii],lineLabel=r'$H_{%u,%u}$' %(n_row[ii]+1,m_row[ii]+1))[0:2]
        
        plt.figure(f_mag_phase.number)
        ax_mag,ax_phase=plot_FRF_mag_phase(f_col,H_cols[:,ii],ax_mag,ax_phase,islin_y=False,lineLabel=r'$H_{%u,%u}$' %(n_row[ii]+1,m_row[ii]+1))[0:2]
    
    f_3d.gca().legend(loc='upper left')
    f_3d.savefig('H_3D.pdf')
    f_Nyq.gca().legend()
    f_Nyq.savefig('H_Nyq.pdf', bbox_inches='tight')
    ax_r.legend()
    f_r_i.savefig('H_r_i.pdf', bbox_inches='tight')
    ax_mag.legend()
    f_mag_phase.savefig('H_mag_phase.pdf', bbox_inches='tight')
    
    plt.show()