function stat_n = occupation_statistics_1d_bhm2C_Fcn(psic,basis)


%%
stat_n = struct();


%%
n_bs = basis.n_bs;
L = basis.L;

state_up = basis.state_up;
state_dn = basis.state_dn;


%% probability list
prob_ltC = abs(conj(psic).*psic);


%%
density_up_ltC = sum(state_up .* repmat(prob_ltC,1,L),1);
parity_up_ltC = sum(mod(state_up,2) .* repmat(prob_ltC,1,L),1);

density_dn_ltC = sum(state_dn .* repmat(prob_ltC,1,L),1);
parity_dn_ltC = sum(mod(state_dn,2) .* repmat(prob_ltC,1,L),1);


%%
stat_n.L = L;
stat_n.n_bs = n_bs;

stat_n.prob_ltC = prob_ltC;

stat_n.density_up_ltC = density_up_ltC;
stat_n.parity_up_ltC = parity_up_ltC;

stat_n.density_dn_ltC = density_dn_ltC;
stat_n.parity_dn_ltC = parity_dn_ltC;

