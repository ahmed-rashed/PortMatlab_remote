function [EigVectors_Normalized, EigValues_vec]=MDOF_Eig_Visc(M_mat, C_mat, K_mat, isPropotional)

if nargin<4
    isPropotional=false;
end

N=size(M_mat,1);
if isPropotional || all(all(C_mat==0))    %Undamped or proportional
    [EigVectors_U,EigValues_U_mat]=eig(-K_mat,M_mat);
    EigValues_U_vec=diag(EigValues_U_mat);
    
    M_r_vec=diag(EigVectors_U.'*M_mat*EigVectors_U);
    w_U_r_vec=sqrt(-EigValues_U_vec);

    C_r_vec=diag(EigVectors_U.'*C_mat*EigVectors_U);
    zeta_r_vec=C_r_vec./(2*M_r_vec.*w_U_r_vec);
    w_d_r_vec=w_U_r_vec.*sqrt(1-zeta_r_vec.^2);
    EigValues_vec_temp=-zeta_r_vec.*w_U_r_vec+1i*w_d_r_vec;

    EigValues_vec=nan(2*N,1);
    EigValues_vec(1:2:2*N-1)=EigValues_vec_temp;
    EigValues_vec(2:2:2*N)=conj(EigValues_vec_temp);
    
    EigVectors_Normalized=nan(N,2*N);
    EigVectors_Normalized(:,1:2:2*N-1)=EigVectors_U*diag(1./sqrt(1i*2*w_d_r_vec.*M_r_vec));
    EigVectors_Normalized(:,2:2:2*N)=conj(EigVectors_Normalized(:,1:2:2*N-1));
else    %Non-proportional
    [EigVectors_Normalized,EigValues_vec]=quad_eig(K_mat,C_mat,M_mat);
end

[~,Index]=sort(abs(imag(EigValues_vec(1:2:2*N-1))));
IIndex=nan(2*N,1);
IIndex(1:2:2*N-1)=2*Index-1;
IIndex(2:2:2*N)=2*Index;
EigValues_vec=EigValues_vec(IIndex);
EigVectors_Normalized=EigVectors_Normalized(:,IIndex);