function [ham_elems] = hamiltonian_1d_bhm_elements(basis,BDC)
% [ham_elems] = hamiltonian_1d_bhm_elements(basis,BDC)
%

if nargin < 2
    BDC = 'obc';
elseif nargin > 2
    error('Error! Invalid input parameters!')
end


%%
ham_elems = struct();


%%
ns = basis.ns;
L = basis.L;
% nMax = basis.nMax;
state = basis.state;


%% diag_term
i_lt = (1:ns)';
j_lt = i_lt;

% part 01: on-site interaction term
fprintf('\ngenerating hamiltonian elements of U.\n')
tVal_U = tic;
coeff_diagTerm_U = state.*(state-1)/2;
ham_elems_U = struct();
for kk = 1:L
    field_cur = ['site_',num2str(kk)];
    ham_elems_U.(field_cur) = sparse(i_lt,j_lt,...
        coeff_diagTerm_U(:,kk),ns,ns);
end
ham_elems_U.sumAllSite = sparse(i_lt,j_lt,...
    sum(coeff_diagTerm_U,2),ns,ns);
tD_U = toc(tVal_U);
fprintf('elapsed time is %.6f seconds.\n',tD_U)


% part 02: chemical potential, |up>
fprintf('\ngenerating hamiltonian elements of mu.\n')
tVal_mu = tic;
ham_elems_mu = struct();
for kk = 1:L
    field_cur = ['site_',num2str(kk)];
    ham_elems_mu.(field_cur) = sparse(i_lt,j_lt,...
        state(:,kk),ns,ns);
end
ham_elems_mu.sumAllSite = sparse(i_lt,j_lt,...
    sum(state,2),ns,ns);
tD_mu = toc(tVal_mu);
fprintf('elapsed time is %.6f seconds.\n',tD_mu)


%% off-diag_trems, tunneling terms
fprintf('\ngenerating hamiltonian elements of J.\n')
t_stJ = tic;
[ham_elems_J] = ...
    tunneling_matrix_elements_1dbhm(basis,BDC);
tDJ = toc(t_stJ);
fprintf('elapsed time is %.6f seconds.\n',tDJ)

   
%% results return
ham_elems.ham_elems_U = ham_elems_U;
ham_elems.ham_elems_J = ham_elems_J;
ham_elems.ham_elems_mu = ham_elems_mu;

ham_elems.BDC = BDC;



