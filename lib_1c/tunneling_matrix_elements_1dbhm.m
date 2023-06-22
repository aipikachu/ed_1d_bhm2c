function [elems_J] = tunneling_matrix_elements_1dbhm(basis,BDC)


%%
elems_J = struct();


%%
ns = basis.ns;
L = basis.L;
nMax = basis.nMax;
state = basis.state;


%% index generate
idx_i_all = struct();
idx_j_all = struct();
idx_k_all = struct();

idx_ilt = [];
idx_jlt = [];
idx_klt = [];

switch BDC
    case 'obc'
        for kk = 1:L-1
            ia = kk;
            ib = ia + 1;
            
            field_cur = ['site_',num2str(ia),'_',num2str(ib)];
            idx_i_all.(field_cur) = [];
            idx_j_all.(field_cur) = [];
            idx_k_all.(field_cur) = [];

        end
        
    case 'pbc'
        for kk = 1:L
            ia = kk;
            ib = mod(kk,L)+1;
            
            field_cur = ['site_',num2str(ia),'_',num2str(ib)];
            idx_i_all.(field_cur) = [];
            idx_j_all.(field_cur) = [];
            idx_k_all.(field_cur) = [];
            
        end
        
    otherwise
        error('Invalid boundary condition!')

end


%%
% notes
% Creation and annihilation operators: Adag, A
seq_all = (1:ns)';
idxstatel = basis.idxIntl;

switch BDC

    case 'obc'
        for jj = 1:L-1
            ia = jj;
            ib = ia + 1;
            
            field_cur = ['site_',num2str(ia),'_',num2str(ib)];
            
            % step_01: for |up> state
            % Adag*A
            idx_selC = (state(:,ia) < nMax) & (state(:,ib) > 0);
            seq_stat_cur = seq_all(idx_selC);

            state_cur = state(idx_selC,:);
            state_nxt = state_cur;
            state_nxt(:,ia) = state_nxt(:,ia) + 1;
            state_nxt(:,ib) = state_nxt(:,ib) - 1;
            idx_lt_nxt = state_nxt * ((nMax+1).^(L-1:-1:0))';           
            [~,seq_stat_nxt] = ismember(idx_lt_nxt,...
                idxstatel);
            coeff_stat = -1 * sqrt(state_cur(:,ia)+1) ...
                .* sqrt(state_cur(:,ib));
            
            %
            idx_i_all.(field_cur) = cat(1,...
                idx_i_all.(field_cur),seq_stat_nxt);
            idx_j_all.(field_cur) = cat(1,...
                idx_j_all.(field_cur),seq_stat_cur);
            idx_k_all.(field_cur) = cat(1,...
                idx_k_all.(field_cur),coeff_stat);
            
            idx_ilt = cat(1,idx_ilt,seq_stat_nxt);
            idx_jlt = cat(1,idx_jlt,seq_stat_cur);
            idx_klt = cat(1,idx_klt,coeff_stat);
            
            
            % A*Adag
            idx_selC = (state(:,ia) > 0) & (state(:,ib) < nMax);
            seq_stat_cur = seq_all(idx_selC);

            state_cur = state(idx_selC,:);
            state_nxt = state_cur;
            state_nxt(:,ia) = state_nxt(:,ia) - 1;
            state_nxt(:,ib) = state_nxt(:,ib) + 1;
            idx_lt_nxt = state_nxt * ((nMax+1).^(L-1:-1:0))';
            [~,seq_stat_nxt] = ismember(idx_lt_nxt,...
                idxstatel);
            coeff_stat = -1 * sqrt(state_cur(:,ia)) ...
                .* sqrt(state_cur(:,ib)+1);
            
            %
            idx_i_all.(field_cur) = cat(1,...
                idx_i_all.(field_cur),seq_stat_nxt);
            idx_j_all.(field_cur) = cat(1,...
                idx_j_all.(field_cur),seq_stat_cur);
            idx_k_all.(field_cur) = cat(1,...
                idx_k_all.(field_cur),coeff_stat);
            
            idx_ilt = cat(1,idx_ilt,seq_stat_nxt);
            idx_jlt = cat(1,idx_jlt,seq_stat_cur);
            idx_klt = cat(1,idx_klt,coeff_stat);
            
        end
        
    case 'pbc'
        for jj = 1:L
            ia = jj;
            ib = mod(jj,L)+1;
            
            field_cur = ['site_',num2str(ia),'_',num2str(ib)];
            
            % step_01: for |up> state
            % Adag*A  
            idx_selC = (state(:,ia) < nMax) & (state(:,ib) > 0);
            seq_stat_cur = seq_all(idx_selC);

            state_cur = state(idx_selC,:);
            state_nxt = state_cur;
            state_nxt(:,ia) = state_nxt(:,ia) + 1;
            state_nxt(:,ib) = state_nxt(:,ib) - 1;
            idx_lt_nxt = state_nxt * ((nMax+1).^(L-1:-1:0))';
            [~,seq_stat_nxt] = ismember(idx_lt_nxt,...
                idxstatel);
            coeff_stat = -1 * sqrt(state_cur(:,ia)+1) ...
                .* sqrt(state_cur(:,ib));
            
            %
            idx_i_all.(field_cur) = cat(1,...
                idx_i_all.(field_cur),seq_stat_nxt);
            idx_j_all.(field_cur) = cat(1,...
                idx_j_all.(field_cur),seq_stat_cur);
            idx_k_all.(field_cur) = cat(1,...
                idx_k_all.(field_cur),coeff_stat);
            
            idx_ilt = cat(1,idx_ilt,seq_stat_nxt);
            idx_jlt = cat(1,idx_jlt,seq_stat_cur);
            idx_klt = cat(1,idx_klt,coeff_stat);
            
            
            % A*Adag
            idx_selC = (state(:,ia) > 0) & (state(:,ib) < nMax);
            seq_stat_cur = seq_all(idx_selC);

            state_cur = state(idx_selC,:);
            state_nxt = state_cur;
            state_nxt(:,ia) = state_nxt(:,ia) - 1;
            state_nxt(:,ib) = state_nxt(:,ib) + 1;
            idx_lt_nxt = state_nxt * ((nMax+1).^(L-1:-1:0))';
            [~,seq_stat_nxt] = ismember(idx_lt_nxt,...
                idxstatel);
            coeff_stat = -1 * sqrt(state_cur(:,ia)) ...
                .* sqrt(state_cur(:,ib)+1);
            
            %
            idx_i_all.(field_cur) = cat(1,...
                idx_i_all.(field_cur),seq_stat_nxt);
            idx_j_all.(field_cur) = cat(1,...
                idx_j_all.(field_cur),seq_stat_cur);
            idx_k_all.(field_cur) = cat(1,...
                idx_k_all.(field_cur),coeff_stat);
            
            idx_ilt = cat(1,idx_ilt,seq_stat_nxt);
            idx_jlt = cat(1,idx_jlt,seq_stat_cur);
            idx_klt = cat(1,idx_klt,coeff_stat);
            

        end  
        
end


%%
switch BDC
    case 'obc'
        for kk = 1:L-1
            ia = kk;
            ib = ia + 1;
            field_cur = ['site_',num2str(ia),'_',num2str(ib)];
            
            i_cur = idx_i_all.(field_cur);
            j_cur = idx_j_all.(field_cur);
            v_cur = idx_k_all.(field_cur);

            elems_J.(field_cur) = sparse(i_cur,j_cur,v_cur,ns,ns);
            
        end
        
    case 'pbc'
        for kk = 1:L
            ia = kk;
            ib = mod(kk,L)+1;
            field_cur = ['site_',num2str(ia),'_',num2str(ib)];
            
            i_cur = idx_i_all.(field_cur);
            j_cur = idx_j_all.(field_cur);
            v_cur = idx_k_all.(field_cur);
            
            elems_J.(field_cur) = sparse(i_cur,j_cur,v_cur,ns,ns);
           
        end  
        
end
    
%
elems_J.sumAllSite = sparse(idx_ilt,idx_jlt,...
    idx_klt,ns,ns);

