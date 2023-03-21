clearvars
clc
close all

f_final=40;
n_f_points=10000;
islin=0;

N=3;
m=10*ones(1,N);
c=8*ones(1,N+1);
k=[161000/2*ones(1,N),0];

ii_row=[1,1,1];
jj_row=[3,2,1];

[M,C,K]=N_DOF_sys(m,c,k);isProportional=false;
% C=zeros(size(M));isProportional=true;
% C=.0001*M+.0001*K;isProportional=true;

f_col=linspace(0,f_final,n_f_points).';
w_col=2*pi*f_col;

%% Slow FRF
tic
H_cols0=MDOF_FRF_Visc_slow(M,C,K,w_col,ii_row,jj_row);
toc

tic
H_cols1=MDOF_FRF_Visc_Slow_Fortran(M,C,K,w_col,ii_row,jj_row);clear MDOF_FRF_Visc_Slow_Fortran;
toc
max(max(abs(H_cols1-H_cols0)))

tic
H_cols2=MDOF_FRF_Visc_Slow_cpp(M,C,K,w_col,ii_row,jj_row);clear MDOF_FRF_Visc_Slow_cpp;
toc
max(max(abs(H_cols2-H_cols0)))

%% Fast FRF
tic
[EigVectors_Normalized0, EigValues_vec0]=MDOF_Eig_Visc(M, C, K,isProportional);
H_cols3=MDOF_FRF_Visc(EigValues_vec0, EigVectors_Normalized0, w_col, ii_row, jj_row);
toc
max(max(abs(H_cols3-H_cols0)))
H_cols=H_cols3;

tic
[EigVectors_Normalized1, EigValues_vec1]=MDOF_Eig_Visc_Fortran(M, C, K,isProportional);clear MDOF_Eig_Visc_Fortran;
H_cols4=MDOF_FRF_Visc_Fortran(EigValues_vec1, EigVectors_Normalized1, w_col, ii_row, jj_row);clear MDOF_FRF_Visc_Fortran;
toc
max(max(abs(EigVectors_Normalized1-EigVectors_Normalized0)))
max(abs(EigValues_vec1-EigValues_vec0))
max(max(abs(H_cols4-H_cols3)))
H_cols=H_cols4;

tic
[EigVectors_Normalized2, EigValues_vec2]=MDOF_Eig_Visc_cpp(M, C, K,isProportional);clear MDOF_Eig_Visc_cpp;
H_cols5=MDOF_FRF_Visc_cpp(EigValues_vec2, EigVectors_Normalized2, w_col, ii_row, jj_row);clear MDOF_FRF_Visc_cpp;
toc
max(max(abs(EigVectors_Normalized2-EigVectors_Normalized0)))
max(abs(EigValues_vec2-EigValues_vec0))
max(max(abs(H_cols5-H_cols3)))
H_cols=H_cols5;

tic
[EigVectors_Normalized3, EigValues_vec3]=MDOF_Eig_Visc_cpp_Armadillo(M, C, K,isProportional);clear MDOF_Eig_Visc_cpp_Armadillo;
H_cols6=MDOF_FRF_Visc_cpp_Armadillo(EigValues_vec3, EigVectors_Normalized3, w_col, ii_row, jj_row);clear MDOF_FRF_Visc_cpp_Armadillo;
toc
max(max(abs(EigVectors_Normalized3-EigVectors_Normalized0)))
max(abs(EigValues_vec3-EigValues_vec0))
max(max(abs(H_cols6-H_cols3)))
H_cols=H_cols6;

f_3d=figure;hold all

f_r_i=figure;
ax_r=subplot(2,1,1);hold all
ax_i=subplot(2,1,2);hold all

f_Nyq=figure;hold all

f_mag_phase=figure;
ax_mag=subplot(4,1,[1,2,3]);hold all
ax_phase=subplot(4,1,4);hold all

legend_str=cell(length(jj_row),1);
for ii=1:length(jj_row)
    legend_str(ii)=cellstr("$H_{1,"+jj_row(ii)+'}$');
    
    figure(f_3d)
    plot_FRF_3d(f_col,H_cols(:,ii),'','',1);

    plot_FRF_r_i(f_col,H_cols(:,ii),ax_r,ax_i);

    figure(f_Nyq)
    plot_FRF_Nyq(H_cols(:,ii));

    plot_FRF_mag_phase(f_col,H_cols(:,ii),islin,ax_mag,ax_phase);
end

figure(f_3d)
legend(legend_str,'interpreter','latex','Location','bestOutside')

figure(f_Nyq)
legend(legend_str,'interpreter','latex','Location','bestOutside')

legend(ax_r,legend_str,'interpreter','latex','Location','NorthEast')

legend(ax_mag,legend_str,'interpreter','latex','Location','NorthEast')

export_figure((1:4),'',["H_MDOF_3D","H_MDOF_r_i","H_MDOF_Nyq","H_MDOF_mag_phase"])