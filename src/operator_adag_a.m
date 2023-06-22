function [op] = operator_adag_a(basis,ia,ib)
%
% return 'Adag_ia*A_ib'

%%
ns = basis.ns;
L = basis.L;
nMax = basis.nMax;
state = basis.state;
idxstatel = basis.idxIntl;


%%
seq_all = (1:ns)';

% Adag*A
idx_selC = (state(:,ia) < nMax) & (state(:,ib) > 0);
j_lt = seq_all(idx_selC);

state_cur = state(idx_selC,:);
state_nxt = state_cur;
state_nxt(:,ia) = state_nxt(:,ia) + 1;
state_nxt(:,ib) = state_nxt(:,ib) - 1;
idx_lt_nxt = state_nxt * ((nMax+1).^(L-1:-1:0))';           
[~,i_lt] = ismember(idx_lt_nxt,idxstatel);
k_lt = 1 * sqrt(state_cur(:,ia)+1) .* sqrt(state_cur(:,ib));


%%
op = sparse(i_lt,j_lt,k_lt,ns,ns);
