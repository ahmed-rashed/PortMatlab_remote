import matplotlib.pyplot as plt

def plot_FRF_r_i(f_vec,H_vec,ax_r=None,ax_i=None,f_label='$f$ (Hz)',H_subtitle='H',lineLabel='',DispMagLines=False):
    import numpy as np

    if ax_r is None:
        ax_r=plt.subplot(2,1,1)
        if ax_i is None:
            plt.setp(ax_r.get_xticklabels(), visible=False)
            
    if ax_i is None:
        ax_i=plt.subplot(2,1,2, sharex=ax_r)

    n_f=f_vec.size
    np.reshape(f_vec,n_f)
    np.reshape(H_vec,n_f)

    ind1=f_label.find('$')
    if ind1==-1:
        raise NameError('f_label does not include LaTeX inline equation !!')
    ind2=f_label.find('$',ind1+1)
    if ind2==-1:
        raise NameError('f_label does not include LaTeX inline equation !!')
    H_subtitle1=r'%s\left(%s\right)' %(H_subtitle ,f_label[ind1+1:ind2])

    h1=ax_r.plot(f_vec,H_vec.real,label=r'%s' %lineLabel)
    ax_r.set_ylabel(r'$\Re\left(%s\right)$' %H_subtitle1)
    ax_r.grid()
    
    h2=ax_i.plot(f_vec,H_vec.imag)
    ax_i.set_xlabel(r'%s' %f_label)
    ax_i.set_ylabel(r'$\Im\left(%s\right)$' %H_subtitle1)
    ax_i.grid()
    
    return ax_r,ax_i,h1,h2