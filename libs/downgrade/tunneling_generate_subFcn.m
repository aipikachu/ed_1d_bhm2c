function [J_elems_info] = tunneling_generate_subFcn(basis,BDC,n_st,n_ed)

%%
J_elems_info = struct();


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
            field_cur = ['site_',num2str(kk),'_',num2str(kk+1)];
            idx_up_i_all.(field_cur) = [];
            idx_up_j_all.(field_cur) = [];
            idx_up_k_all.(field_cur) = [];
            
            idx_dn_i_all.(field_cur) = [];
            idx_dn_j_all.(field_cur) = [];
            idx_dn_k_all.(field_cur) = [];

        end
        
    case 'pbc'
        for kk = 1:L
            field_cur = ['site_',num2str(kk),'_',num2str(mod(kk,L)+1)];
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


%% off-diag_trems, tunneling terms
for kk = n_st:n_ed
    if mod(kk,100) == 1
        fprintf('current process: %d / %d, ~ %.6f%s.\n',...
            kk,n_bs,(kk/n_bs)*100,'%')
    end
    
    stat_up_curr = state_up(kk,:);
    stat_dn_curr = state_dn(kk,:);
    
    stat_idx_cur = kk;
    
    % notes
    % |up> : A, Adag
    % |dn> : B, Bdag 
    
    switch BDC

        case 'obc'
            for jj = 1:L-1
                field_cur = ['site_',num2str(jj),'_',num2str(jj+1)];
                
                % step_01: for |up> state
                % Adag*A
                if stat_up_curr(jj+1)
                    stat_up_nxt = stat_up_curr;
                    stat_up_nxt(jj) = stat_up_nxt(jj) + 1;
                    stat_up_nxt(jj+1) = stat_up_nxt(jj+1) - 1;
                    
                    if ~sum(stat_up_nxt>nMax)
                        stat_idx_nxt = state_index_search(stat_up_nxt,stat_dn_curr,basis);
                        stat_up_coeff_nxt = -1 * sqrt(stat_up_curr(jj)+1) ...
                            * sqrt(stat_up_curr(jj+1));
                        
                        %                    
                        idx_up_i_all.(field_cur) = cat(1,...
                            idx_up_i_all.(field_cur),stat_idx_nxt);
                        idx_up_j_all.(field_cur) = cat(1,...
                            idx_up_j_all.(field_cur),stat_idx_cur);
                        idx_up_k_all.(field_cur) = cat(1,...
                            idx_up_k_all.(field_cur),stat_up_coeff_nxt);
                        
                        idx_up_ilt = cat(1,idx_up_ilt,stat_idx_nxt);
                        idx_up_jlt = cat(1,idx_up_jlt,stat_idx_cur);
                        idx_up_klt = cat(1,idx_up_klt,stat_up_coeff_nxt);
                        
                    end
                end

                % A*Adag
                if stat_up_curr(jj)
                    stat_up_nxt = stat_up_curr;
                    stat_up_nxt(jj) = stat_up_nxt(jj) - 1;
                    stat_up_nxt(jj+1) = stat_up_nxt(jj+1) + 1;
                    
                    if ~sum(stat_up_nxt>nMax)
                        stat_idx_nxt = state_index_search(stat_up_nxt,stat_dn_curr,basis);
                        stat_up_coeff_nxt = -1 * sqrt(stat_up_curr(jj)) ...
                            * sqrt(stat_up_curr(jj+1)+1);
                        
                        % 
                        idx_up_i_all.(field_cur) = cat(1,...
                            idx_up_i_all.(field_cur),stat_idx_nxt);
                        idx_up_j_all.(field_cur) = cat(1,...
                            idx_up_j_all.(field_cur),stat_idx_cur);
                        idx_up_k_all.(field_cur) = cat(1,...
                            idx_up_k_all.(field_cur),stat_up_coeff_nxt);
                        
                        idx_up_ilt = cat(1,idx_up_ilt,stat_idx_nxt);
                        idx_up_jlt = cat(1,idx_up_jlt,stat_idx_cur);
                        idx_up_klt = cat(1,idx_up_klt,stat_up_coeff_nxt);
                        
                    end
                end
                
                % step_02: for |dn> state
                % Bdag*B
                if stat_dn_curr(jj+1)
                    stat_dn_nxt = stat_dn_curr;
                    stat_dn_nxt(jj) = stat_dn_nxt(jj) + 1;
                    stat_dn_nxt(jj+1) = stat_dn_nxt(jj+1) - 1;
                    
                    if ~sum(stat_dn_nxt>nMax)
                        stat_idx_nxt = state_index_search(stat_up_curr,stat_dn_nxt,basis);
                        stat_dn_coeff_nxt = -1 * sqrt(stat_dn_curr(jj)+1) ...
                            * sqrt(stat_dn_curr(jj+1));
                        
                        % 
                        idx_dn_i_all.(field_cur) = cat(1,...
                            idx_dn_i_all.(field_cur),stat_idx_nxt);
                        idx_dn_j_all.(field_cur) = cat(1,...
                            idx_dn_j_all.(field_cur),stat_idx_cur);
                        idx_dn_k_all.(field_cur) = cat(1,...
                            idx_dn_k_all.(field_cur),stat_dn_coeff_nxt);
                        
                        idx_dn_ilt = cat(1,idx_dn_ilt,stat_idx_nxt);
                        idx_dn_jlt = cat(1,idx_dn_jlt,stat_idx_cur);
                        idx_dn_klt = cat(1,idx_dn_klt,stat_dn_coeff_nxt);
                        
                    end
                end

                % B*Bdag
                if stat_dn_curr(jj)
                    stat_dn_nxt = stat_dn_curr;
                    stat_dn_nxt(jj) = stat_dn_nxt(jj) - 1;
                    stat_dn_nxt(jj+1) = stat_dn_nxt(jj+1) + 1;
                    
                    if ~sum(stat_dn_nxt>nMax)
                        stat_idx_nxt = state_index_search(stat_up_curr,stat_dn_nxt,basis);
                        stat_dn_coeff_nxt = -1 * sqrt(stat_dn_curr(jj)) ...
                            * sqrt(stat_dn_curr(jj+1)+1);
                        
                        % 
                        idx_dn_i_all.(field_cur) = cat(1,...
                            idx_dn_i_all.(field_cur),stat_idx_nxt);
                        idx_dn_j_all.(field_cur) = cat(1,...
                            idx_dn_j_all.(field_cur),stat_idx_cur);
                        idx_dn_k_all.(field_cur) = cat(1,...
                            idx_dn_k_all.(field_cur),stat_dn_coeff_nxt);
                        
                        idx_dn_ilt = cat(1,idx_dn_ilt,stat_idx_nxt);
                        idx_dn_jlt = cat(1,idx_dn_jlt,stat_idx_cur);
                        idx_dn_klt = cat(1,idx_dn_klt,stat_dn_coeff_nxt);
                        
                    end
                end
                
            end
            
        case 'pbc'
            for jj = 1:L
                field_cur = ['site_',num2str(jj),'_',num2str(mod(jj,L)+1)];
                
                % step_01: for |up> state
                % Adag*A
                if stat_up_curr(mod(jj,L)+1)
                    stat_up_nxt = stat_up_curr;
                    stat_up_nxt(jj) = stat_up_nxt(jj) + 1;
                    stat_up_nxt(mod(jj,L)+1) = stat_up_nxt(mod(jj,L)+1) - 1;
                    
                    if ~sum(stat_up_nxt>nMax)
                        stat_idx_nxt = state_index_search(stat_up_nxt,stat_dn_curr,basis);
                        stat_up_coeff_nxt = -1 * sqrt(stat_up_curr(jj)+1) ...
                            * sqrt(stat_up_curr(mod(jj,L)+1));
                        
                        %
                        idx_up_i_all.(field_cur) = cat(1,...
                            idx_up_i_all.(field_cur),stat_idx_nxt);
                        idx_up_j_all.(field_cur) = cat(1,...
                            idx_up_j_all.(field_cur),stat_idx_cur);
                        idx_up_k_all.(field_cur) = cat(1,...
                            idx_up_k_all.(field_cur),stat_up_coeff_nxt);
                        
                        idx_up_ilt = cat(1,idx_up_ilt,stat_idx_nxt);
                        idx_up_jlt = cat(1,idx_up_jlt,stat_idx_cur);
                        idx_up_klt = cat(1,idx_up_klt,stat_up_coeff_nxt);
                        
                    end
                end

                % A*Adag
                if stat_up_curr(jj)
                    stat_up_nxt = stat_up_curr;
                    stat_up_nxt(jj) = stat_up_nxt(jj) - 1;
                    stat_up_nxt(mod(jj,L)+1) = stat_up_nxt(mod(jj,L)+1) + 1;
                    
                    if ~sum(stat_up_nxt>nMax)
                        stat_idx_nxt = state_index_search(stat_up_nxt,stat_dn_curr,basis);
                        stat_up_coeff_nxt = -1 * sqrt(stat_up_curr(jj)) ...
                            * sqrt(stat_up_curr(mod(jj,L)+1)+1);
                        
                        % 
                        idx_up_i_all.(field_cur) = cat(1,...
                            idx_up_i_all.(field_cur),stat_idx_nxt);
                        idx_up_j_all.(field_cur) = cat(1,...
                            idx_up_j_all.(field_cur),stat_idx_cur);
                        idx_up_k_all.(field_cur) = cat(1,...
                            idx_up_k_all.(field_cur),stat_up_coeff_nxt);
                        
                        idx_up_ilt = cat(1,idx_up_ilt,stat_idx_nxt);
                        idx_up_jlt = cat(1,idx_up_jlt,stat_idx_cur);
                        idx_up_klt = cat(1,idx_up_klt,stat_up_coeff_nxt);
                        
                    end
                end
                
                % step_02: for |dn> state
                % Bdag*B
                if stat_dn_curr(mod(jj,L)+1)
                    stat_dn_nxt = stat_dn_curr;
                    stat_dn_nxt(jj) = stat_dn_nxt(jj) + 1;
                    stat_dn_nxt(mod(jj,L)+1) = stat_dn_nxt(mod(jj,L)+1) - 1;
                    
                    if ~sum(stat_dn_nxt>nMax)
                        stat_idx_nxt = state_index_search(stat_up_curr,stat_dn_nxt,basis);
                        stat_dn_coeff_nxt = -1 * sqrt(stat_dn_curr(jj)+1) ...
                            * sqrt(stat_dn_curr(mod(jj,L)+1));
                        
                        % 
                        idx_dn_i_all.(field_cur) = cat(1,...
                            idx_dn_i_all.(field_cur),stat_idx_nxt);
                        idx_dn_j_all.(field_cur) = cat(1,...
                            idx_dn_j_all.(field_cur),stat_idx_cur);
                        idx_dn_k_all.(field_cur) = cat(1,...
                            idx_dn_k_all.(field_cur),stat_dn_coeff_nxt);
                        
                        idx_dn_ilt = cat(1,idx_dn_ilt,stat_idx_nxt);
                        idx_dn_jlt = cat(1,idx_dn_jlt,stat_idx_cur);
                        idx_dn_klt = cat(1,idx_dn_klt,stat_dn_coeff_nxt);
                        
                    end
                end

                % B*Bdag
                if stat_dn_curr(jj)
                    stat_dn_nxt = stat_dn_curr;
                    stat_dn_nxt(jj) = stat_dn_nxt(jj) - 1;
                    stat_dn_nxt(mod(jj,L)+1) = stat_dn_nxt(mod(jj,L)+1) + 1;
                    
                    if ~sum(stat_dn_nxt>nMax)
                        stat_idx_nxt = state_index_search(stat_up_curr,stat_dn_nxt,basis);
                        stat_dn_coeff_nxt = -1 * sqrt(stat_dn_curr(jj)) ...
                            * sqrt(stat_dn_curr(mod(jj,L)+1)+1);
                        
                        % 
                        idx_dn_i_all.(field_cur) = cat(1,...
                            idx_dn_i_all.(field_cur),stat_idx_nxt);
                        idx_dn_j_all.(field_cur) = cat(1,...
                            idx_dn_j_all.(field_cur),stat_idx_cur);
                        idx_dn_k_all.(field_cur) = cat(1,...
                            idx_dn_k_all.(field_cur),stat_dn_coeff_nxt);
                        
                        idx_dn_ilt = cat(1,idx_dn_ilt,stat_idx_nxt);
                        idx_dn_jlt = cat(1,idx_dn_jlt,stat_idx_cur);
                        idx_dn_klt = cat(1,idx_dn_klt,stat_dn_coeff_nxt);
                        
                    end
                end
                
            end
      
    end

end


%%
J_elems_info.idx_up_i_all = idx_up_i_all;
J_elems_info.idx_up_j_all = idx_up_j_all;
J_elems_info.idx_up_k_all = idx_up_k_all;
J_elems_info.idx_dn_i_all = idx_dn_i_all;
J_elems_info.idx_dn_j_all = idx_dn_j_all;
J_elems_info.idx_dn_k_all = idx_dn_k_all;
J_elems_info.idx_up_ilt = idx_up_ilt;
J_elems_info.idx_up_jlt = idx_up_jlt;
J_elems_info.idx_up_klt = idx_up_klt;
J_elems_info.idx_dn_ilt = idx_dn_ilt;
J_elems_info.idx_dn_jlt = idx_dn_jlt;
J_elems_info.idx_dn_klt = idx_dn_klt;

