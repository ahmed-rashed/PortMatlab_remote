function [Epsi_normalized,Val_vec]=quad_eig(K_mat,C_mat,M_mat)

N=size(M_mat,1);

N_mat=M_mat;
A=[-K_mat,0*M_mat;0*M_mat,N_mat];
B=[C_mat,M_mat;N_mat,0*M_mat];

[Phi,Val_mat] = eig(A,B);
Val_vec=diag(Val_mat);

B_r_vec=diag(Phi.'*B*Phi);

Phi_normalized=Phi;
for ii=1:2*N
    Phi_normalized(:,ii)=Phi_normalized(:,ii)/sqrt(B_r_vec(ii));
end

Epsi_normalized=Phi_normalized(1:N,:);