function [op] = operator_n(basis,ia)
%
% return n_ia

%
%%
ns = basis.ns;
% L = basis.L;
% nMax = basis.nMax;
state = basis.state;

%%
i_lt = (1:ns)';
j_lt = i_lt;
k_lt = state(:,ia);

%%
op = sparse(i_lt,j_lt,k_lt,ns,ns);
