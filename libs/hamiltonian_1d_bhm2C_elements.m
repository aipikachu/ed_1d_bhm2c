function [ham_elems] = hamiltonian_1d_bhm2C_elements(basis,BDC)
% [ham_elems] = hamiltonian_1d_bhm2C_elements(basis,BDC)
%

if nargin < 2
    BDC = 'obc';
end

%%
ham_elems = struct();


%% Input parameter check
TAG_s = true;
try
    cache_dir = [pwd,filesep,'cache'];
    cache_info = load([cache_dir,filesep,'ham_elems_cache.mat']);
    if (~strcmp(cache_info.BDC,BDC)) ...
            || (cache_info.basis.nMax ~= basis.nMax) ...
            || (cache_info.basis.L ~= basis.L) ...
            || (cache_info.basis.N_up ~= basis.N_up) ...
            || (cache_info.basis.N_dn ~= basis.N_dn)
        TAG_s = false; 
    end
catch
    TAG_s = false;    
end

%
if TAG_s
    ham_elems = cache_info.ham_elems;
    fprintf('Using the cache data!\n')
    return
else
    fprintf('Re-generate the hamiltonian elements!\n')
end


%%
n_bs = basis.n_bs;
L = basis.L;
nMax = basis.nMax;

state_up = basis.state_up;
state_dn = basis.state_dn;

% N_up = basis.N_up;
% N_dn = basis.N_dn;
% NA = N_up + N_dn;


%% diag_term
% part 01_1: on-site interaction term of (up,up)
% coeff_diagTerm_U_uu = sum(state_up.*(state_up-1)/2,2);
% ham_elems_U_uu = sparse(diag(coeff_diagTerm_U_uu));
coeff_diagTerm_U_uu = state_up.*(state_up-1)/2;
ham_elems_U_uu = struct();
for kk = 1:L
    field_cur = ['site_',num2str(kk)];
    ham_elems_U_uu.(field_cur) = sparse(diag(coeff_diagTerm_U_uu(:,kk)));
end
ham_elems_U_uu.sumAllSite = sparse(diag(sum(coeff_diagTerm_U_uu,2)));


% part 01_2: on-site interaction term of (dn,dn)
% coeff_diagTerm_U_dd = sum(state_dn.*(state_dn-1)/2,2);
% ham_elems_U_dd = sparse(diag(coeff_diagTerm_U_dd));
coeff_diagTerm_U_dd = state_dn.*(state_dn-1)/2;
ham_elems_U_dd = struct();
for kk = 1:L
    field_cur = ['site_',num2str(kk)];
    ham_elems_U_dd.(field_cur) = sparse(diag(coeff_diagTerm_U_dd(:,kk)));
end
ham_elems_U_dd.sumAllSite = sparse(diag(sum(coeff_diagTerm_U_dd,2)));


% part 01_3: on-site interaction term of (up,dn)
% coeff_diagTerm_U_ud = sum(state_up.*state_dn,2);
% ham_elems_U_ud = sparse(diag(coeff_diagTerm_U_ud));
coeff_diagTerm_U_ud = state_up.*state_dn;
ham_elems_U_ud = struct();
for kk = 1:L
    field_cur = ['site_',num2str(kk)];
    ham_elems_U_ud.(field_cur) = sparse(diag(coeff_diagTerm_U_ud(:,kk)));
end
ham_elems_U_ud.sumAllSite = sparse(diag(sum(coeff_diagTerm_U_ud,2)));


% part 02_1: chemical potential, |up>
ham_elems_mu_up = struct();
for kk = 1:L
    field_cur = ['site_',num2str(kk)];
    ham_elems_mu_up.(field_cur) = sparse(diag(state_up(:,kk)));
end
ham_elems_mu_up.sumAllSite = sparse(diag(sum(state_up,2)));


% part 02_2: chemical potential, |dn>
ham_elems_mu_dn = struct();
for kk = 1:L
    field_cur = ['site_',num2str(kk)];
    ham_elems_mu_dn.(field_cur) = sparse(diag(state_dn(:,kk)));
end
ham_elems_mu_dn.sumAllSite = sparse(diag(sum(state_dn,2)));


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
% idxlt = basis.idxstatel;

for kk = 1:n_bs
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


%
ham_elems_J_up = struct();
ham_elems_J_dn = struct();

switch BDC
    case 'obc'
        for kk = 1:L-1
            field_cur = ['site_',num2str(kk),'_',num2str(kk+1)];
            
            i_up_cur = idx_up_i_all.(field_cur);
            j_up_cur = idx_up_j_all.(field_cur);
            v_up_cur = idx_up_k_all.(field_cur);
            
            i_dn_cur = idx_dn_i_all.(field_cur);
            j_dn_cur = idx_dn_j_all.(field_cur);
            v_dn_cur = idx_dn_k_all.(field_cur);
            
            ham_elems_J_up.(field_cur) = ...
                sparse(i_up_cur,j_up_cur,v_up_cur,...
                n_bs,n_bs);
            ham_elems_J_dn.(field_cur) = ...
                sparse(i_dn_cur,j_dn_cur,v_dn_cur,...
                n_bs,n_bs);
            
        end
        
    case 'pbc'
        for kk = 1:L
            field_cur = ['site_',num2str(kk),'_',num2str(mod(kk,L)+1)];
            
            i_up_cur = idx_up_i_all.(field_cur);
            j_up_cur = idx_up_j_all.(field_cur);
            v_up_cur = idx_up_k_all.(field_cur);
            
            i_dn_cur = idx_dn_i_all.(field_cur);
            j_dn_cur = idx_dn_j_all.(field_cur);
            v_dn_cur = idx_dn_k_all.(field_cur);
            
            ham_elems_J_up.(field_cur) = ...
                sparse(i_up_cur,j_up_cur,v_up_cur,...
                n_bs,n_bs);
            ham_elems_J_dn.(field_cur) = ...
                sparse(i_dn_cur,j_dn_cur,v_dn_cur,...
                n_bs,n_bs);
            
        end  
        
end
    
%
ham_elems_J_up.sumAllSite = sparse(idx_up_ilt,idx_up_jlt,...
    idx_up_klt,n_bs,n_bs);
ham_elems_J_dn.sumAllSite = sparse(idx_dn_ilt,idx_dn_jlt,...
    idx_dn_klt,n_bs,n_bs);


%% results return
ham_elems.ham_elems_U_uu = ham_elems_U_uu;
ham_elems.ham_elems_U_dd = ham_elems_U_dd;
ham_elems.ham_elems_U_ud = ham_elems_U_ud;
ham_elems.ham_elems_J_up = ham_elems_J_up;
ham_elems.ham_elems_J_dn = ham_elems_J_dn;
ham_elems.ham_elems_mu_up = ham_elems_mu_up;
ham_elems.ham_elems_mu_dn = ham_elems_mu_dn;

ham_elems.BDC = BDC;


%% write cache info. into disk
try
    cache_dir = [pwd,filesep,'cache'];
    if ~exist(cache_dir,'dir')
        mkdir(cache_dir)
    end
    
    save([cache_dir,filesep,'ham_elems_cache.mat'],...
        'basis','BDC','ham_elems');
    
catch
    
end

