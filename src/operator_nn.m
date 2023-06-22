function [op] = operator_nn(basis,ia,ib)
%
% return 'n_ia * n_ib

%
%%
ns = basis.ns;
% L = basis.L;
% nMax = basis.nMax;
state = basis.state;

%%
i_lt = (1:ns)';
j_lt = i_lt;
k_lt = state(:,ia).*state(:,ib);

%%
op = sparse(i_lt,j_lt,k_lt,ns,ns);
