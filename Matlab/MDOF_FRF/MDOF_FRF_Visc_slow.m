function H_cols=MDOF_FRF_Visc_slow(M_mat,C_mat,K_mat,w_col,ii_row,jj_row)

n_f_points=size(w_col,1);
H_cols=zeros(n_f_points,length(jj_row));
for ii=1:n_f_points
    H_w=inv(-w_col(ii)^2*M_mat+1i*w_col(ii)*C_mat+K_mat);

    ind_row=sub2ind(size(M_mat),ii_row,jj_row);
    H_cols(ii,:)=H_w(ind_row);
end