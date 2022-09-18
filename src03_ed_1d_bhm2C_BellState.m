% a example code for generating Bell pairs

close all
clc

addpath(genpath('libs'))


%% init state
L = 6;
st_init_up = [1 0 1 0 1 0];
st_init_dn = [0 1 0 1 0 1];

if (length(st_init_up) ~= L) || (length(st_init_dn) ~= L)
    error('Error! Invalid input parameters!')
end

%
N_up = sum(st_init_up);
N_dn = sum(st_init_dn);
% the maximum occupation number of up/dn state in each lattice site
nMax = 3; 


%% Hubbard parameters
J = 50.0 * 2 * pi;
U = 1000.0 * 2 * pi;

%
U_uu = U * ones(1,L);      % U_up_up
U_dd = U * ones(1,L);      % U_dn_dn
U_ud = U * ones(1,L);      % U_up_dn
J_up = J * mod(1:L-1,2);
J_dn = J * mod(1:L-1,2);

% staggered potential
staG0 = 0.0 * U;  

mu_up_lt = 0.0*2*pi*(1:L) + mod(1:L,2)*staG0;
mu_dn_lt = 0.0*2*pi*(1:L) + mod(1:L,2)*staG0;

BDC = 'obc';   % 'obc' - open boundary condition; 'pbs' - periodic case
% BDC = 'pbc';   % 'obc' - open boundary condition; 'pbs' - periodic case


%% evolve time for quantum walk
J_ex = 4*(J/2/pi)^2/(U/2/pi);
nt = 101;
% T = 10 * 1e-03;          % 10 ms
T = 1/J_ex/4;          % 10 ms
tl = linspace(0,T,nt);
dt = tl(2) - tl(1);


%% basis generate
basis = boson_basis_spinor(L,N_up,N_dn,nMax);
ns = basis.n_bs;
fprintf('Total basis number is %d.\n',ns)


%% state index search
idx_init_st = state_index_search(st_init_up,st_init_dn,basis);

% basis.state_up(idx_init_st,:)
% basis.state_dn(idx_init_st,:)

phi_init = zeros(ns,1);
phi_init(idx_init_st) = 1;


%% hamiltonian elements generate
tStart = tic; 
ham_elems = hamiltonian_1d_bhm2C_elements(basis,BDC);
tEnd = toc(tStart);
fprintf('Elapsed time of hamiltonian elements generation is %.6f seconds.\n',tEnd)

ham_cur = hamiltonian_1d_bhm2C(basis,ham_elems,J_up,J_dn,...
    U_uu,U_dd,U_ud,mu_up_lt,mu_dn_lt);

ham = ham_cur.ham;
fprintf('Number of nonzero hamiltonian elements: %d.\n',nnz(ham))

% return
%%
density_up_Mt = [];
density_dn_Mt = [];

probl_Mt = [];

psic = phi_init;
probl = abs(conj(psic).*psic);
probl_Mt = cat(2,probl_Mt,probl);

stat_nC = occupation_statistics_1d_bhm2C_Fcn(psic,basis);
density_up_Mt = cat(1,density_up_Mt,stat_nC.density_up_ltC);
density_dn_Mt = cat(1,density_dn_Mt,stat_nC.density_dn_ltC);

%
tStart = tic; 
for kk = 2:nt
    fprintf('Current process: %04d / %04d.\n',kk,nt)
    
    tSC = tic;
    psic = expv(-1i*dt,ham,psic,1e-8,50);
    tEC = toc(tSC);
    fprintf('time for evolution: %.6f seconds.\n',tEC)
    
    stat_nC = occupation_statistics_1d_bhm2C_Fcn(psic,basis);
    density_up_Mt = cat(1,density_up_Mt,stat_nC.density_up_ltC);
    density_dn_Mt = cat(1,density_dn_Mt,stat_nC.density_dn_ltC);
    
    probl = abs(conj(psic).*psic);
    probl_Mt = cat(2,probl_Mt,probl);
    
    % fprintf('\n')
end
tEnd = toc(tStart);
fprintf('\nElapsed time is %.6f seconds.\n',tEnd)


%%
x = 1:L;
y = tl*1000;

% density profile of the |up> state
figure('Color','w','Position',[120 120 560 420])
imagesc(x,y,density_up_Mt)
ax = gca;
ax.FontSize = 14;
xlabel('Sites','FontSize',16)
ylabel('Evolution time (ms)','FontSize',16)


% density profile of the |dn> state
figure('Color','w','Position',[720 120 560 420])
imagesc(x,y,density_dn_Mt)
ax = gca;
ax.FontSize = 14;
xlabel('Sites','FontSize',16)
ylabel('Evolution time (ms)','FontSize',16)


%% probablity of final state
psi_fs = psic;
probl_fs = abs(conj(psi_fs).*psi_fs);


figure('Color','w','Position',[420 220 560 420])
plot(probl_fs)
xlim([1,ns])
ax = gca;
ax.FontSize = 14;
xlabel('Basis number','FontSize',16)
ylabel('Probability','FontSize',16)


%
[probl_fs_sort,idx_sort] = sort(probl_fs,'descend');
k_s = 8;
fprintf('\nShow the first %d largest state as follows:\n',k_s)
for kk = 1:k_s
    stat_up_cur = basis.state_up(idx_sort(kk),:);
    stat_dn_cur = basis.state_dn(idx_sort(kk),:);
    
    fprintf('%04d) The %d-th state with p_c = %.6f\n',...
        kk,kk,probl_fs_sort(kk))
    fprintf('\t The current |up> state: %s\n',num2str(stat_up_cur))
    fprintf('\t The current |dn> state: %s\n',num2str(stat_dn_cur))
    
end
fprintf('\n')


%%
% figure('Color','w','Position',[720 120 560 420])
% imagesc(probl_Mt)
% ax = gca;
% ax.FontSize = 14;
% xlabel('Sites','FontSize',16)
% ylabel('Evolution time (ms)','FontSize',16)
% % probl_Mt


