function [basis] = boson_basis_spinor(L,N_up,N_dn,nMax)

%%
% L = 4;
% N_up = 1;
% N_dn = 3;
% nMax = 3;

%%
basis = struct();

%%
if nargin < 4
    nMax = 2;
end

basis_up = boson_basis_1d(L,N_up,nMax);
basis_dn = boson_basis_1d(L,N_dn,nMax);

n_up = basis_up.ns;
n_dn = basis_dn.ns;

n_bs = n_up * n_dn;

% for state up
state_up = reshape(repmat(reshape(basis_up.state,1,[]),n_dn,1),[],L);
idxIntl_up = reshape(repmat(reshape(basis_up.idxIntl,1,[]),n_dn,1),[],1);

% for state down
state_dn = repmat(basis_dn.state,n_up,1);
idxIntl_dn = repmat(basis_dn.idxIntl,n_up,1);

% idx list
idxstatel = [idxIntl_up,idxIntl_dn];


%%
basis.state_up = state_up;
basis.state_dn = state_dn;
basis.n_bs = n_bs;
basis.idxstatel = idxstatel;
basis.idxIntl_up = idxIntl_up;
basis.idxIntl_dn = idxIntl_dn;
basis.n_up = n_up;
basis.n_dn = n_dn;

basis.nMax = nMax;
basis.L = L;
basis.N_up = N_up;
basis.N_dn = N_dn;

