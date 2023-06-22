function [s] = obs_1dbhm_N(psi,basis)
%% [] = obs_1dbhm_N(psi_in,basis)
% state representation
% |psi> = \sum_i \alphi_i |\alphi_i>


%%
s = struct();


%%
L = basis.L;
% ns = basis.ns;
% nMax = basis.nMax;
% state = basis.state;


%% observation of 'density' and 'parity' distribution
% tic
for kk = 1:L
    op_nc = operator_n(basis,kk);

    s.density.site_n(kk) = kk;
    s.density.vals(kk) = psi'*op_nc*psi;

    s.parity.site_n(kk) = kk;
    s.parity.vals(kk) = psi'*mod(op_nc,2)*psi;

end
% toc

