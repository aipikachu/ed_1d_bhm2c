function [s] = obs_1dbhm_NN(psi,basis)
%% [] = obs_1dbhm_NN(psi_in,basis)
% state representation: 
% |psi> = \sum_i \alphi_i |\alphi_i>


%%
s = struct();


%%
L = basis.L;
% ns = basis.ns;
% nMax = basis.nMax;
% state = basis.state;


%% initial
tic
% 1) <ni*nj>
s.obs_ninj.nn_site = zeros(L*L,2);
s.obs_ninj.nn_vals = zeros(L*L,1);
% s.obs_ninj.dist_label  = zeros()
% 2) (<ni*nj>-<ni><nj>)
s.obs_corr.nn_site = zeros(L*L,2);
s.obs_corr.nn_vals = zeros(L*L,1);
toc


%%
nn_dist = struct();
co_dist = struct();
for kk = 1:2*L-1
    nn_dist.(['dist_',num2str(kk)]).sites = [];
    nn_dist.(['dist_',num2str(kk)]).vals = [];
    nn_dist.(['dist_',num2str(kk)]).num_d = 0;

    co_dist.(['dist_',num2str(kk)]).sites = [];
    co_dist.(['dist_',num2str(kk)]).vals = [];
    co_dist.(['dist_',num2str(kk)]).num_d = 0;
end


%%
tic
for ia = 1:L
    for ib = 1:L
        kk = ib + (ia-1)*L;

        dist_cur = ib - ia;
        dist_abs = dist_cur + L;

        % operators
        op_ninj = operator_nn(basis,ia,ib);

        val_ninj = psi'*op_ninj*psi;
        val_corr = abs(val_ninj ...
            - (psi'*operator_n(basis,ia)*psi) ...
            * (psi'*operator_n(basis,ib)*psi));

        % for <ni*nj>
        s.obs_ninj.nn_site(kk,:) = [ia,ib];
        s.obs_ninj.nn_vals(kk) = val_ninj;

        % for (<ni*nj>-<ni><nj>)
        s.obs_corr.nn_site(kk,:) = [ia,ib];
        s.obs_corr.nn_vals(kk) = val_corr;

        %
        nn_dist.(['dist_',num2str(dist_abs)]).sites = ...
            cat(1,nn_dist.(['dist_',num2str(dist_abs)]).sites,[ia,ib]);
        nn_dist.(['dist_',num2str(dist_abs)]).vals = ...
            cat(1,nn_dist.(['dist_',num2str(dist_abs)]).vals,val_ninj);
        nn_dist.(['dist_',num2str(dist_abs)]).num_d = ...
            nn_dist.(['dist_',num2str(dist_abs)]).num_d + 1;

        %
        co_dist.(['dist_',num2str(dist_abs)]).sites = ...
            cat(1,co_dist.(['dist_',num2str(dist_abs)]).sites,[ia,ib]);
        co_dist.(['dist_',num2str(dist_abs)]).vals = ...
            cat(1,co_dist.(['dist_',num2str(dist_abs)]).vals,val_corr);
        co_dist.(['dist_',num2str(dist_abs)]).num_d = ...
            co_dist.(['dist_',num2str(dist_abs)]).num_d + 1;

    end
end
toc

s.obs_ninj.nn_dist = nn_dist;
s.obs_corr.co_dist = co_dist;


%%
nn_vs_dist = zeros(1,2*L-1);
co_vs_dist = zeros(1,2*L-1);
nn_sum = 0;
co_sum = 0;
for kk = 1:2*L-1
    nn_vs_dist(kk) = sum(nn_dist.(['dist_',num2str(kk)]).vals) ...
        / nn_dist.(['dist_',num2str(kk)]).num_d;
    co_vs_dist(kk) = sum(co_dist.(['dist_',num2str(kk)]).vals) ...
        / co_dist.(['dist_',num2str(kk)]).num_d;

    nn_sum = nn_sum + sum(nn_dist.(['dist_',num2str(kk)]).vals);
    co_sum = co_sum + sum(co_dist.(['dist_',num2str(kk)]).vals);

end

s.obs_ninj.nn_vs_dist = nn_vs_dist;
s.obs_ninj.dist_list = 1-L:L-1;
s.obs_ninj.nn_avg_val = nn_sum/L/L;
s.obs_corr.co_vs_dist = co_vs_dist;
s.obs_corr.dist_list = 1-L:L-1;
s.obs_corr.co_avg_val = co_sum/L/L;


return
%%
figure('Color','w')
hold on
plot(1:2*L-1,co_vs_dist)
plot(1:2*L-1,nn_vs_dist)
hold off
