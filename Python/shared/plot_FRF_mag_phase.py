from correctedPhase import correctedPhase

def plot_FRF_mag_phase(f_vec,H_vec,ax_mag=None,ax_phase=None,maxPhaseLag=None,islin_y=True,islin_x=True,f_label='$f$ (Hz)',H_subtitle='H',lineLabel='',DispMagLines=False):
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator
    from matplotlib import rc
    
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)

    if ax_mag is None:
        ax_mag=plt.subplot2grid((4,1), (0,0), rowspan=3)
        if ax_phase is None:
            plt.setp(ax_mag.get_xticklabels(), visible=False)
            
    if ax_phase is None:
        ax_phase=plt.subplot2grid((4,1), (3,0), sharex=ax_mag)

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

    h1=ax_mag.plot(f_vec,np.absolute(H_vec),label=r'%s' %lineLabel)
    
    if DispMagLines:
        n_MagLines=20
        delta_temp=np.floor(n_f/n_MagLines)
        iidx=np.arange(1,n_f,delta_temp,dtype=np.intp)
    
        markerline, stemlines, baseline=ax_mag.stem(f_vec[iidx],np.absolute(H_vec[iidx]),linestyle=line_3d[0].get_linestyle())
        plt.setp(baseline,'linewidth', h1[0].get_linewidth()/3)
        ax_mag.plot(f_vec[iidx],np.absolute(H_vec[iidx]),'.',markersize=15)
    
    if islin_y:
        ax_mag.set_yscale('linear')
    else:
        ax_mag.set_yscale('log')

    if islin_x:
        ax_mag.set_xscale('linear')
        #ax_phase.set_xscale('linear')
    else:
        ax_mag.set_xscale('log')
        #ax_phase.set_xscale('log')
    
    ax_mag.set_ylabel(r'$\left|%s\right|$' %H_subtitle1)
    ax_mag.grid()
    
    H_angle_vec=correctedPhase(H_vec,maxPhaseLag)
    h2=ax_phase.plot(f_vec,H_angle_vec)
    if DispMagLines:
        ax_phase.plot(f_vec[iidx],H_angle_vec[iidx],'.',markersize=15)
    
    ax_phase.set_xlabel(r'%s' %f_label)
    ax_phase.set_ylabel(r'$\angle %s$' %H_subtitle1)
       
    yStep=np.pi
    if ax_phase.get_autoscaley_on():
        yLim_min=yStep*np.floor(np.amin(H_angle_vec)/yStep)
        yLim_max=yStep*np.ceil(np.amax(H_angle_vec)/yStep)
    else:
        ylimits=ax_phase.get_ylim()
        yLim_min=np.minimum(yStep*np.floor(np.amin(H_angle_vec)/yStep),ylimits[0])
        yLim_max=np.maximum(yStep*np.ceil(np.amax(H_angle_vec)/yStep),ylimits[1])
    
    N_ticks=np.round((yLim_max-yLim_min)/yStep)+1
    if N_ticks>6:
        yStep=yStep*np.round(N_ticks/6)

    yTicks=np.arange(yLim_min,yLim_max+yStep,yStep)
    N_ticks=yTicks.size
    ax_phase.set_yticks(yTicks)

    
    yTickLabels=[r'%d\pi' %(yTicks[ii]/np.pi) for ii in np.arange(0,N_ticks)]
    for ii in np.arange(0,N_ticks):
        if yTicks[ii]==0:
            yTickLabels[ii]='0'
        elif yTicks[ii]==np.pi:
            yTickLabels[ii]=r'\pi'
        elif yTicks[ii]==-np.pi:
            yTickLabels[ii]=r'-\pi'

    ax_phase.set_yticklabels(yTickLabels)
    if yLim_min!=yLim_max:
        ax_phase.set_ylim(yLim_min,yLim_max)

    ax_phase.grid()
    if yStep!=np.pi:
        ax_phase.yaxis.set_minor_locator(MultipleLocator(np.pi))
        ax_phase.tick_params(axis='y',which='minor')
    
    return ax_mag,ax_phase,h1,h2