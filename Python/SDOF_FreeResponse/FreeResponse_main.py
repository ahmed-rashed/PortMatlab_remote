from numpy import linspace,sqrt,pi,array
import matplotlib.pyplot as plt
import matplotlib
from Free_Response import Free_Response

if __name__ == '__main__':
    plt.close('all')
    matplotlib.rcParams['text.usetex'] = True
    
    x0=-1
    v0=0
    w_n=1
    f_n=w_n/2.0/pi
    T_n=1/f_n
    
    zeta_vec=array([0,.1,.2,.4,1/sqrt(2),1,2])
    
    t_vec=linspace(0,2*T_n,500)
    
    for zeta in zeta_vec:
        x_vec=Free_Response(w_n,zeta,x0,v0,t_vec)
        plt.plot(t_vec/T_n,x_vec, label=r'$\zeta = %.3g$' %zeta)
    
    plt.title(r'$x(t)$ for $\omega_{n}=%g$, $x_{0}=%g$ and $v_{0}=%g$' %(w_n,x0,v0))
    plt.xlabel(r'$t/T_{n}\qquad,:T_{n}=1/f_{n}=2\pi/\omega_{n}$')
    
    plt.legend(loc='lower right')
    plt.grid()
    
    plt.savefig('SDOF_1_Python.pdf', bbox_inches='tight')
    
    plt.show()