function [elems_Jup,elems_Jdn] = tunneling_matrix_elements_generate(basis,BDC)


%%
elems_Jup = struct();
elems_Jdn = struct();


%%
n_bs = basis.n_bs;
L = basis.L;
nMax = basis.nMax;
state_up = basis.state_up;
state_dn = basis.state_dn;


%% index generate
idx_up_i_all = struct();
idx_up_j_all = struct();
idx_up_k_all = struct();
idx_dn_i_all = struct();
idx_dn_j_all = struct();
idx_dn_k_all = struct();
idx_up_ilt = [];
idx_up_jlt = [];
idx_up_klt = [];
idx_dn_ilt = [];
idx_dn_jlt = [];
idx_dn_klt = [];

switch BDC
    case 'obc'
        for kk = 1:L-1
            ia = kk;
            ib = ia + 1;
            
            field_cur = ['site_',num2str(ia),'_',num2str(ib)];
            idx_up_i_all.(field_cur) = [];
            idx_up_j_all.(field_cur) = [];
            idx_up_k_all.(field_cur) = [];
            
            idx_dn_i_all.(field_cur) = [];
            idx_dn_j_all.(field_cur) = [];
            idx_dn_k_all.(field_cur) = [];

        end
        
    case 'pbc'
        for kk = 1:L
            ia = kk;
            ib = mod(kk,L)+1;
            
            field_cur = ['site_',num2str(ia),'_',num2str(ib)];
            idx_up_i_all.(field_cur) = [];
            idx_up_j_all.(field_cur) = [];
            idx_up_k_all.(field_cur) = [];
            
            idx_dn_i_all.(field_cur) = [];
            idx_dn_j_all.(field_cur) = [];
            idx_dn_k_all.(field_cur) = [];
            
        end
        
    otherwise
        error('Invalid boundary condition!')

end


%%
% notes
% |up> : A, Adag
% |dn> : B, Bdag 
seq_all = (1:n_bs)';
idxstatel = basis.idxstatel;
idxlt_up = basis.idxIntl_up;
idxlt_dn = basis.idxIntl_dn;

switch BDC

    case 'obc'
        for jj = 1:L-1
            ia = jj;
            ib = ia + 1;
            
            field_cur = ['site_',num2str(ia),'_',num2str(ib)];
            
            % step_01: for |up> state
            % Adag*A
            idx_selC = (state_up(:,ia) < nMax) & (state_up(:,ib) > 0);
            seq_stat_cur = seq_all(idx_selC);
            
            idx_lt_dn_nxt = idxlt_dn(idx_selC);
            state_up_cur = state_up(idx_selC,:);
            state_up_nxt = state_up_cur;
            state_up_nxt(:,ia) = state_up_nxt(:,ia) + 1;
            state_up_nxt(:,ib) = state_up_nxt(:,ib) - 1;
            idx_lt_up_nxt = state_up_nxt * ((nMax+1).^(L-1:-1:0))';           
            [~,seq_stat_nxt] = ismember([idx_lt_up_nxt,idx_lt_dn_nxt],...
                idxstatel,'rows');
            coeff_up_stat = -1 * sqrt(state_up_cur(:,ia)+1) ...
                .* sqrt(state_up_cur(:,ib));
            
            %
            idx_up_i_all.(field_cur) = cat(1,...
                idx_up_i_all.(field_cur),seq_stat_nxt);
            idx_up_j_all.(field_cur) = cat(1,...
                idx_up_j_all.(field_cur),seq_stat_cur);
            idx_up_k_all.(field_cur) = cat(1,...
                idx_up_k_all.(field_cur),coeff_up_stat);
            
            idx_up_ilt = cat(1,idx_up_ilt,seq_stat_nxt);
            idx_up_jlt = cat(1,idx_up_jlt,seq_stat_cur);
            idx_up_klt = cat(1,idx_up_klt,coeff_up_stat);
            
            
            % A*Adag
            idx_selC = (state_up(:,ia) > 0) & (state_up(:,ib) < nMax);
            seq_stat_cur = seq_all(idx_selC);
            
            idx_lt_dn_nxt = idxlt_dn(idx_selC);
            state_up_cur = state_up(idx_selC,:);
            state_up_nxt = state_up_cur;
            state_up_nxt(:,ia) = state_up_nxt(:,ia) - 1;
            state_up_nxt(:,ib) = state_up_nxt(:,ib) + 1;
            idx_lt_up_nxt = state_up_nxt * ((nMax+1).^(L-1:-1:0))';
            [~,seq_stat_nxt] = ismember([idx_lt_up_nxt,idx_lt_dn_nxt],...
                idxstatel,'rows');
            coeff_up_stat = -1 * sqrt(state_up_cur(:,ia)) ...
                .* sqrt(state_up_cur(:,ib)+1);
            
            %
            idx_up_i_all.(field_cur) = cat(1,...
                idx_up_i_all.(field_cur),seq_stat_nxt);
            idx_up_j_all.(field_cur) = cat(1,...
                idx_up_j_all.(field_cur),seq_stat_cur);
            idx_up_k_all.(field_cur) = cat(1,...
                idx_up_k_all.(field_cur),coeff_up_stat);
            
            idx_up_ilt = cat(1,idx_up_ilt,seq_stat_nxt);
            idx_up_jlt = cat(1,idx_up_jlt,seq_stat_cur);
            idx_up_klt = cat(1,idx_up_klt,coeff_up_stat);
            

            % step_02: for |dn> state
            % Bdag*B
            idx_selC = (state_dn(:,ia) < nMax) & (state_dn(:,ib) > 0);
            seq_stat_cur = seq_all(idx_selC);
            
            idx_lt_up_nxt = idxlt_up(idx_selC);
            state_dn_cur = state_dn(idx_selC,:);
            state_dn_nxt = state_dn_cur;
            state_dn_nxt(:,ia) = state_dn_nxt(:,ia) + 1;
            state_dn_nxt(:,ib) = state_dn_nxt(:,ib) - 1;
            idx_lt_dn_nxt = state_dn_nxt * ((nMax+1).^(L-1:-1:0))';
            [~,seq_stat_nxt] = ismember([idx_lt_up_nxt,idx_lt_dn_nxt],...
                idxstatel,'rows');
            coeff_dn_stat = -1 * sqrt(state_dn_cur(:,ia)+1) ...
                .* sqrt(state_dn_cur(:,ib));
            
            %
            idx_dn_i_all.(field_cur) = cat(1,...
                idx_dn_i_all.(field_cur),seq_stat_nxt);
            idx_dn_j_all.(field_cur) = cat(1,...
                idx_dn_j_all.(field_cur),seq_stat_cur);
            idx_dn_k_all.(field_cur) = cat(1,...
                idx_dn_k_all.(field_cur),coeff_dn_stat);
            
            idx_dn_ilt = cat(1,idx_dn_ilt,seq_stat_nxt);
            idx_dn_jlt = cat(1,idx_dn_jlt,seq_stat_cur);
            idx_dn_klt = cat(1,idx_dn_klt,coeff_dn_stat);
            
            
            % B*Bdag
            idx_selC = (state_dn(:,ia) > 0) & (state_dn(:,ib) < nMax);
            seq_stat_cur = seq_all(idx_selC);
            
            idx_lt_up_nxt = idxlt_up(idx_selC);
            state_dn_cur = state_dn(idx_selC,:);
            state_dn_nxt = state_dn_cur;
            state_dn_nxt(:,ia) = state_dn_nxt(:,ia) - 1;
            state_dn_nxt(:,ib) = state_dn_nxt(:,ib) + 1;
            idx_lt_dn_nxt = state_dn_nxt * ((nMax+1).^(L-1:-1:0))';
            [~,seq_stat_nxt] = ismember([idx_lt_up_nxt,idx_lt_dn_nxt],...
                idxstatel,'rows');
            coeff_dn_stat = -1 * sqrt(state_dn_cur(:,ia)) ...
                .* sqrt(state_dn_cur(:,ib)+1);
            
            %
            idx_dn_i_all.(field_cur) = cat(1,...
                idx_dn_i_all.(field_cur),seq_stat_nxt);
            idx_dn_j_all.(field_cur) = cat(1,...
                idx_dn_j_all.(field_cur),seq_stat_cur);
            idx_dn_k_all.(field_cur) = cat(1,...
                idx_dn_k_all.(field_cur),coeff_dn_stat);
            
            idx_dn_ilt = cat(1,idx_dn_ilt,seq_stat_nxt);
            idx_dn_jlt = cat(1,idx_dn_jlt,seq_stat_cur);
            idx_dn_klt = cat(1,idx_dn_klt,coeff_dn_stat);

        end
        
    case 'pbc'
        for jj = 1:L
            ia = jj;
            ib = mod(jj,L)+1;
            
            field_cur = ['site_',num2str(ia),'_',num2str(ib)];
            
            % step_01: for |up> state
            % Adag*A  
            idx_selC = (state_up(:,ia) < nMax) & (state_up(:,ib) > 0);
            seq_stat_cur = seq_all(idx_selC);
            
            idx_lt_dn_nxt = idxlt_dn(idx_selC);
            state_up_cur = state_up(idx_selC,:);
            state_up_nxt = state_up_cur;
            state_up_nxt(:,ia) = state_up_nxt(:,ia) + 1;
            state_up_nxt(:,ib) = state_up_nxt(:,ib) - 1;
            idx_lt_up_nxt = state_up_nxt * ((nMax+1).^(L-1:-1:0))';
            [~,seq_stat_nxt] = ismember([idx_lt_up_nxt,idx_lt_dn_nxt],...
                idxstatel,'rows');
            coeff_up_stat = -1 * sqrt(state_up_cur(:,ia)+1) ...
                .* sqrt(state_up_cur(:,ib));
            
            %
            idx_up_i_all.(field_cur) = cat(1,...
                idx_up_i_all.(field_cur),seq_stat_nxt);
            idx_up_j_all.(field_cur) = cat(1,...
                idx_up_j_all.(field_cur),seq_stat_cur);
            idx_up_k_all.(field_cur) = cat(1,...
                idx_up_k_all.(field_cur),coeff_up_stat);
            
            idx_up_ilt = cat(1,idx_up_ilt,seq_stat_nxt);
            idx_up_jlt = cat(1,idx_up_jlt,seq_stat_cur);
            idx_up_klt = cat(1,idx_up_klt,coeff_up_stat);
            
            
            % A*Adag
            idx_selC = (state_up(:,ia) > 0) & (state_up(:,ib) < nMax);
            seq_stat_cur = seq_all(idx_selC);
            
            idx_lt_dn_nxt = idxlt_dn(idx_selC);
            state_up_cur = state_up(idx_selC,:);
            state_up_nxt = state_up_cur;
            state_up_nxt(:,ia) = state_up_nxt(:,ia) - 1;
            state_up_nxt(:,ib) = state_up_nxt(:,ib) + 1;
            idx_lt_up_nxt = state_up_nxt * ((nMax+1).^(L-1:-1:0))';
            [~,seq_stat_nxt] = ismember([idx_lt_up_nxt,idx_lt_dn_nxt],...
                idxstatel,'rows');
            coeff_up_stat = -1 * sqrt(state_up_cur(:,ia)) ...
                .* sqrt(state_up_cur(:,ib)+1);
            
            %
            idx_up_i_all.(field_cur) = cat(1,...
                idx_up_i_all.(field_cur),seq_stat_nxt);
            idx_up_j_all.(field_cur) = cat(1,...
                idx_up_j_all.(field_cur),seq_stat_cur);
            idx_up_k_all.(field_cur) = cat(1,...
                idx_up_k_all.(field_cur),coeff_up_stat);
            
            idx_up_ilt = cat(1,idx_up_ilt,seq_stat_nxt);
            idx_up_jlt = cat(1,idx_up_jlt,seq_stat_cur);
            idx_up_klt = cat(1,idx_up_klt,coeff_up_stat);
            

            % step_02: for |dn> state
            % Bdag*B
            idx_selC = (state_dn(:,ia) < nMax) & (state_dn(:,ib) > 0);
            seq_stat_cur = seq_all(idx_selC);
            
            idx_lt_up_nxt = idxlt_up(idx_selC);
            state_dn_cur = state_dn(idx_selC,:);
            state_dn_nxt = state_dn_cur;
            state_dn_nxt(:,ia) = state_dn_nxt(:,ia) + 1;
            state_dn_nxt(:,ib) = state_dn_nxt(:,ib) - 1;
            idx_lt_dn_nxt = state_dn_nxt * ((nMax+1).^(L-1:-1:0))';
            [~,seq_stat_nxt] = ismember([idx_lt_up_nxt,idx_lt_dn_nxt],...
                idxstatel,'rows');
            coeff_dn_stat = -1 * sqrt(state_dn_cur(:,ia)+1) ...
                .* sqrt(state_dn_cur(:,ib));
            
            %
            idx_dn_i_all.(field_cur) = cat(1,...
                idx_dn_i_all.(field_cur),seq_stat_nxt);
            idx_dn_j_all.(field_cur) = cat(1,...
                idx_dn_j_all.(field_cur),seq_stat_cur);
            idx_dn_k_all.(field_cur) = cat(1,...
                idx_dn_k_all.(field_cur),coeff_dn_stat);
            
            idx_dn_ilt = cat(1,idx_dn_ilt,seq_stat_nxt);
            idx_dn_jlt = cat(1,idx_dn_jlt,seq_stat_cur);
            idx_dn_klt = cat(1,idx_dn_klt,coeff_dn_stat);
            
            
            % B*Bdag
            idx_selC = (state_dn(:,ia) > 0) & (state_dn(:,ib) < nMax);
            seq_stat_cur = seq_all(idx_selC);
            
            idx_lt_up_nxt = idxlt_up(idx_selC);
            state_dn_cur = state_dn(idx_selC,:);
            state_dn_nxt = state_dn_cur;
            state_dn_nxt(:,ia) = state_dn_nxt(:,ia) - 1;
            state_dn_nxt(:,ib) = state_dn_nxt(:,ib) + 1;
            idx_lt_dn_nxt = state_dn_nxt * ((nMax+1).^(L-1:-1:0))';
            [~,seq_stat_nxt] = ismember([idx_lt_up_nxt,idx_lt_dn_nxt],...
                idxstatel,'rows');
            coeff_dn_stat = -1 * sqrt(state_dn_cur(:,ia)) ...
                .* sqrt(state_dn_cur(:,ib)+1);
            
            %
            idx_dn_i_all.(field_cur) = cat(1,...
                idx_dn_i_all.(field_cur),seq_stat_nxt);
            idx_dn_j_all.(field_cur) = cat(1,...
                idx_dn_j_all.(field_cur),seq_stat_cur);
            idx_dn_k_all.(field_cur) = cat(1,...
                idx_dn_k_all.(field_cur),coeff_dn_stat);
            
            idx_dn_ilt = cat(1,idx_dn_ilt,seq_stat_nxt);
            idx_dn_jlt = cat(1,idx_dn_jlt,seq_stat_cur);
            idx_dn_klt = cat(1,idx_dn_klt,coeff_dn_stat);

        end  
        
end


%%
switch BDC
    case 'obc'
        for kk = 1:L-1
            ia = kk;
            ib = ia + 1;
            field_cur = ['site_',num2str(ia),'_',num2str(ib)];
            
            i_up_cur = idx_up_i_all.(field_cur);
            j_up_cur = idx_up_j_all.(field_cur);
            v_up_cur = idx_up_k_all.(field_cur);
            
            i_dn_cur = idx_dn_i_all.(field_cur);
            j_dn_cur = idx_dn_j_all.(field_cur);
            v_dn_cur = idx_dn_k_all.(field_cur);
            
            elems_Jup.(field_cur) = ...
                sparse(i_up_cur,j_up_cur,v_up_cur,...
                n_bs,n_bs);
            elems_Jdn.(field_cur) = ...
                sparse(i_dn_cur,j_dn_cur,v_dn_cur,...
                n_bs,n_bs);
            
        end
        
    case 'pbc'
        for kk = 1:L
            ia = kk;
            ib = mod(kk,L)+1;
            field_cur = ['site_',num2str(ia),'_',num2str(ib)];
            
            i_up_cur = idx_up_i_all.(field_cur);
            j_up_cur = idx_up_j_all.(field_cur);
            v_up_cur = idx_up_k_all.(field_cur);
            
            i_dn_cur = idx_dn_i_all.(field_cur);
            j_dn_cur = idx_dn_j_all.(field_cur);
            v_dn_cur = idx_dn_k_all.(field_cur);
            
            elems_Jup.(field_cur) = ...
                sparse(i_up_cur,j_up_cur,v_up_cur,...
                n_bs,n_bs);
            elems_Jdn.(field_cur) = ...
                sparse(i_dn_cur,j_dn_cur,v_dn_cur,...
                n_bs,n_bs);
            
        end  
        
end
    
%
elems_Jup.sumAllSite = sparse(idx_up_ilt,idx_up_jlt,...
    idx_up_klt,n_bs,n_bs);
elems_Jdn.sumAllSite = sparse(idx_dn_ilt,idx_dn_jlt,...
    idx_dn_klt,n_bs,n_bs);


