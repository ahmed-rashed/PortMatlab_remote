from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def plot_FRF_3d(f_vec,H_vec,f_label='$f$ (Hz)',H_subtitle='H',lineLabel='',DispRealImag=False,DispMagLines=False):
    import numpy as np

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

    line_3d=plt.plot(f_vec,H_vec.real,H_vec.imag,label=r'%s' %lineLabel)
    
#    plt.axis('equal')
    plt.grid()

    plt.xlabel(r'%s' %f_label)
    plt.ylabel(r'$\Re\left(%s\right)$' %H_subtitle1)
    plt.gca().set_zlabel(r'$\Im\left(%s\right)$' %H_subtitle1)
    
    xlimm=plt.xlim()
    ylimm=plt.ylim()
    zlimm=plt.gca().get_zlim()
    #Bug workaround
    a1=xlimm[0]
    b1=xlimm[1]
    a3=zlimm[0]
    b3=zlimm[1]
    
    if DispRealImag:
        plt.plot(H_vec.real,H_vec.imag,zs=xlimm[0],zdir='x', color=line_3d[0].get_color(),linestyle=line_3d[0].get_linestyle(),linewidth=line_3d[0].get_linewidth()/3)
        plt.plot(f_vec,H_vec.real,zs=a3,zdir='z', color=line_3d[0].get_color(),linestyle=line_3d[0].get_linestyle(),linewidth=line_3d[0].get_linewidth()/3)
        plt.plot(f_vec,H_vec.imag,zs=ylimm[1],zdir='y', color=line_3d[0].get_color(),linestyle=line_3d[0].get_linestyle(),linewidth=line_3d[0].get_linewidth()/3)
        #Bug Workaround
        plt.xlim(a1,b1)
        plt.ylim(ylimm)
        plt.gca().set_zlim(a3,b3)
        
        
    if DispMagLines:
        n_MagLines=20
        delta_temp=np.floor(n_f/n_MagLines)
        iidx=np.arange(1,n_f,delta_temp,dtype=np.intp)
        
        plt.quiver3D(f_vec[iidx],0*f_vec[iidx],0*f_vec[iidx],0*f_vec[iidx],H_vec[iidx].real,H_vec[iidx].imag, length=.001)
        
        #plt.plot(f_vec[iidx],H_vec[iidx].real,H_vec[iidx].imag,'.','MarkerSize',15,color=line_3d[0].get_color());
        
    return line_3d