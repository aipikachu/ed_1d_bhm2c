% a example code for single atom quantum walk

close all
clc

addpath(genpath('libs'))


%% init state
L = 9;
st_init_up = zeros(1,L);
st_init_dn = mod(1:L,2);
st_init_dn(5) = 0;

if (length(st_init_up) ~= L) || (length(st_init_dn) ~= L)
    error('Error! Invalid input parameters!')
end

%
N_up = sum(st_init_up);
N_dn = sum(st_init_dn);
% the maximum occupation number of up/dn state in each lattice site
nMax = 3; 


%% Hubbard parameters
J = 55.0 * 2 * pi;
U = 780.0 * 2 * pi;

%
U_uu = U;      % U_up_up
U_dd = U;      % U_dn_dn
U_ud = U;      % U_up_dn
J_up = J;
J_dn = J;

% staggered potential
staG0 = 0.5 * U - 2.0 * J;  

mu_up_lt = 0.0*2*pi*(1:L) + mod(1:L,2)*staG0;
mu_dn_lt = 60.0*2*pi*(1:L) + mod(1:L,2)*staG0;

BDC = 'obc';   % 'obc' - open boundary condition; 'pbs' - periodic case
% BDC = 'pbc';   % 'obc' - open boundary condition; 'pbs' - periodic case


%% evolve time for quantum walk
nt = 101;
T = 100 * 1e-03;          % 100 ms
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


%%
density_up_Mt = [];
density_dn_Mt = [];

psic = phi_init;
% probl = abs(conj(psic).*psic);

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
colormap(hot)
colorbar
ax = gca;
ax.FontSize = 14;
xlabel('Sites','FontSize',16)
ylabel('Evolution time (ms)','FontSize',16)


% density profile of the |dn> state
figure('Color','w','Position',[720 120 560 420])
imagesc(x,y,density_dn_Mt)
colormap(hot)
colorbar
ax = gca;
ax.FontSize = 14;
xlabel('Sites','FontSize',16)
ylabel('Evolution time (ms)','FontSize',16)


%% probablity of final state
psi_fs = psic;
probl_fs = abs(conj(psi_fs).*psi_fs);

