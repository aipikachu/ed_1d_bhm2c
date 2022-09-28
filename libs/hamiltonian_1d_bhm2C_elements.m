function [ham_elems] = hamiltonian_1d_bhm2C_elements(basis,BDC,TAG_)
% [ham_elems] = hamiltonian_1d_bhm2C_elements(basis,BDC)
%

if nargin < 2
    BDC = 'obc';
    TAG_ = '';
elseif nargin < 3
    TAG_ = '';
elseif nargin > 3
    error('Error! Invalid input parameters!')
end


%%
ham_elems = struct();


%% Input parameter check
TAG_s = true;

if isempty(TAG_)
    cache_name = 'ham_elems_cache.mat';
else
    cache_name = ['ham_elems_cache_',TAG_,...
        sprintf('L%02d_Nup%02d_Ndn%02d_nMax%01d',...
        basis.L,basis.N_up,basis.N_dn,basis.nMax),'.mat'];
end

try
    cache_dir = [pwd,filesep,'cache'];
    cache_info = load([cache_dir,filesep,cache_name]);
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
% nMax = basis.nMax;
state_up = basis.state_up;
state_dn = basis.state_dn;


%% diag_term
i_lt = (1:n_bs)';
j_lt = i_lt;

% part 01_1: on-site interaction term of (up,up)
fprintf('\ngenerating hamiltonian elements of U_uu.\n')
tVal_U_uu = tic;
coeff_diagTerm_U_uu = state_up.*(state_up-1)/2;
ham_elems_U_uu = struct();
for kk = 1:L
    field_cur = ['site_',num2str(kk)];
    ham_elems_U_uu.(field_cur) = sparse(i_lt,j_lt,...
        coeff_diagTerm_U_uu(:,kk),n_bs,n_bs);
end
ham_elems_U_uu.sumAllSite = sparse(i_lt,j_lt,...
    sum(coeff_diagTerm_U_uu,2),n_bs,n_bs);
tD_U_uu = toc(tVal_U_uu);
fprintf('elapsed time is %.6f seconds.\n',tD_U_uu)


% part 01_2: on-site interaction term of (dn,dn)
fprintf('\ngenerating hamiltonian elements of U_dd.\n')
tVal_U_dd = tic;
coeff_diagTerm_U_dd = state_dn.*(state_dn-1)/2;
ham_elems_U_dd = struct();
for kk = 1:L
    field_cur = ['site_',num2str(kk)];
    ham_elems_U_dd.(field_cur) = sparse(i_lt,j_lt,...
        coeff_diagTerm_U_dd(:,kk),n_bs,n_bs);
end
ham_elems_U_dd.sumAllSite = sparse(i_lt,j_lt,...
    sum(coeff_diagTerm_U_dd,2),n_bs,n_bs);
tD_U_dd = toc(tVal_U_dd);
fprintf('elapsed time is %.6f seconds.\n',tD_U_dd)


% part 01_3: on-site interaction term of (up,dn)
fprintf('\ngenerating hamiltonian elements of U_ud.\n')
tVal_U_ud = tic;
coeff_diagTerm_U_ud = state_up.*state_dn;
ham_elems_U_ud = struct();
for kk = 1:L
    field_cur = ['site_',num2str(kk)];
    ham_elems_U_ud.(field_cur) = sparse(i_lt,j_lt,...
        coeff_diagTerm_U_ud(:,kk),n_bs,n_bs);
end
ham_elems_U_ud.sumAllSite = sparse(i_lt,j_lt,...
    sum(coeff_diagTerm_U_ud,2),n_bs,n_bs);
tD_U_ud = toc(tVal_U_ud);
fprintf('elapsed time is %.6f seconds.\n',tD_U_ud)


% part 02_1: chemical potential, |up>
fprintf('\ngenerating hamiltonian elements of mu_up.\n')
tVal_mu_up = tic;
ham_elems_mu_up = struct();
for kk = 1:L
    field_cur = ['site_',num2str(kk)];
    ham_elems_mu_up.(field_cur) = sparse(i_lt,j_lt,...
        state_up(:,kk),n_bs,n_bs);
end
ham_elems_mu_up.sumAllSite = sparse(i_lt,j_lt,...
    sum(state_up,2),n_bs,n_bs);
tD_mu_up = toc(tVal_mu_up);
fprintf('elapsed time is %.6f seconds.\n',tD_mu_up)


% part 02_2: chemical potential, |dn>
fprintf('\ngenerating hamiltonian elements of mu_dn.\n')
tVal_mu_dn = tic;
ham_elems_mu_dn = struct();
for kk = 1:L
    field_cur = ['site_',num2str(kk)];
    ham_elems_mu_dn.(field_cur) = sparse(i_lt,j_lt,...
        state_dn(:,kk),n_bs,n_bs);
end
ham_elems_mu_dn.sumAllSite = sparse(i_lt,j_lt,...
    sum(state_dn,2),n_bs,n_bs);
tD_mu_dn = toc(tVal_mu_dn);
fprintf('elapsed time is %.6f seconds.\n',tD_mu_dn)


%% off-diag_trems, tunneling terms
fprintf('\ngenerating hamiltonian elements of J_up and J_dn.\n')
t_stJ = tic;
[ham_elems_J_up,ham_elems_J_dn] = ...
    tunneling_matrix_elements_generate(basis,BDC);
tDJ = toc(t_stJ);
fprintf('elapsed time is %.6f seconds.\n',tDJ)

   
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
    
    save([cache_dir,filesep,cache_name],...
        'basis','BDC','ham_elems','TAG_');
    
catch
    
end

