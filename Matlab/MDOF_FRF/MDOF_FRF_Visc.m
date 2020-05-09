function H_w_n_m_cols=MDOF_FRF_Visc(EigValues_vec, EigVectors_Normalized, w_column, n_row, m_row)

N=size(EigVectors_Normalized,1);
N_w=size(w_column,1);
if (any(size(n_row)~=size(m_row)));error('Dimensions of n_row and m_row must be identical');end
N_cols=size(n_row,2);

H_w_n_m_cols_SDOF=zeros(N_w,N_cols);
H_w_n_m_cols=zeros(N_w,N_cols);
A_ind_row=sub2ind([N,N],n_row,m_row);
for r=1:2*N
    A_r=EigVectors_Normalized(:,r)*EigVectors_Normalized(:,r).';
    A_r_temp_row=A_r(A_ind_row);
    H_w_n_m_cols_SDOF=H_w_n_m_cols_SDOF+(1./(1i*w_column-EigValues_vec(r)))*A_r_temp_row;
    
    if imag(EigValues_vec(r))~=0 && mod(r,2)~=0   %complex eigenvalue and odd r
        continue
    end
    
    H_w_n_m_cols=H_w_n_m_cols+H_w_n_m_cols_SDOF;
    H_w_n_m_cols_SDOF=zeros(N_w,N_cols);
end