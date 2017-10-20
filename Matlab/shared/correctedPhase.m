function H_angle_vec=correctedPhase(H_vec, ...
                                    maxPhaseLag) %Optional arguments
if nargin<2
    maxPhaseLag=-pi;
else
    if isempty(maxPhaseLag)
        maxPhaseLag=-pi;
    end
end

n_f=length(H_vec);

H_angle_vec=angle(H_vec);

H_angle_vec(abs(H_angle_vec)<1e3*eps)=0;
H_angle_vec(abs(H_angle_vec-pi)<1e3*eps)=pi;
H_angle_vec(abs(H_angle_vec+pi)<1e3*eps)=-pi;
H_angle_vec(abs(H_vec)==0)=nan;

%Force phase to be less than maxPhaseLag
indd=H_angle_vec-maxPhaseLag>0;
H_angle_vec(indd)=H_angle_vec(indd)-2*pi*ceil((H_angle_vec(indd)-maxPhaseLag)/2/pi);

%Force phase to be decreasing for jumps >= 45 deg
for ii=1:n_f-1
    dd=H_angle_vec(ii+1)-H_angle_vec(ii);
    if dd>=pi/4
        H_angle_vec(ii+1)=H_angle_vec(ii+1)-2*pi*ceil(dd/2/pi);
    end
end

H_angle_vec=unwrap(H_angle_vec);