import numpy

def correctedPhase(H_vec,maxPhaseLag=-numpy.pi):
    from numpy import pi
    import numpy as np
    
    if maxPhaseLag is None:
        maxPhaseLag=-numpy.pi
        
    n_f=H_vec.size
    np.reshape(H_vec,n_f)
    
    H_angle_vec=np.angle(H_vec)
    
    H_angle_vec[np.absolute(H_angle_vec)<1e3*np.finfo(type(H_vec[0])).eps]=0
    H_angle_vec[np.absolute(H_angle_vec-pi)<1e3*np.finfo(type(H_vec[0])).eps]=pi
    H_angle_vec[np.absolute(H_angle_vec+pi)<1e3*np.finfo(type(H_vec[0])).eps]=-pi
    H_angle_vec[np.absolute(H_vec)==0]=np.nan
    
    #Force phase to be less than maxPhaseLag
    indd=H_angle_vec-maxPhaseLag>0
    H_angle_vec[indd]=H_angle_vec[indd]-2*pi*np.ceil((H_angle_vec[indd]-maxPhaseLag)/2/pi)
    
    #Force phase to be decreasing for jumps >= 45 deg
    for ii in range(0,n_f-1):
        dd=H_angle_vec[ii+1]-H_angle_vec[ii]
        if dd>=pi/4:
            H_angle_vec[ii+1]=H_angle_vec[ii+1]-2*pi*np.ceil(dd/2/pi)
    
    H_angle_vec=np.unwrap(H_angle_vec)
    
    return H_angle_vec