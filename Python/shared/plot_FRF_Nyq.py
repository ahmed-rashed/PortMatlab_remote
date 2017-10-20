import matplotlib.pyplot as plt

def plot_FRF_Nyq(H_vec,H_subtitle='H',lineLabel='',DispMagLines=False):
    import numpy as np

    n_f=H_vec.size
    np.reshape(H_vec,n_f)
    
    H_subtitle1=r'%s\left(f \right)' %H_subtitle
    
    h3=plt.plot(H_vec.real,H_vec.imag,label=r'%s' %lineLabel)

    if DispMagLines:
        n_MagLines=20
        delta_temp=np.floor(n_f/n_MagLines)
        iidx=np.arange(1,n_f,delta_temp,dtype=np.intp)

        plt.quiver(0*H_vec,0*H_vec,H_vec.real,H_vec.imag)
    
    plt.axis('equal')
    
    plt.xlabel(r'$\Re\left(%s\right)$' %H_subtitle1)
    plt.ylabel(r'$\Im\left(%s\right)$' %H_subtitle1)

    plt.grid()
    
    return h3