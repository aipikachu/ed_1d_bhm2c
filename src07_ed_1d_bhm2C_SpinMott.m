% a example code of spin-mott simulation

close all
clc

addpath(genpath('libs'))

%% init state
L = 5;
% st_init_up = [1 1 1 1 1];
% st_init_dn = [1 1 1 1 1];
st_init_up = ones(1,L);
st_init_dn = ones(1,L);

if (length(st_init_up) ~= L) || (length(st_init_dn) ~= L)
    error('Error! Invalid input parameters!')
end

%
N_up = sum(st_init_up);
N_dn = sum(st_init_dn);
% the maximum occupation number of up/dn state in each lattice site
nMax = 3; 


%% Part 01
unit_tilT = 220 * 2 *  pi;
unit_VL = 360 * 2 *  pi;

%
J = 0.085 * unit_tilT;
U = 3.765 * unit_tilT;
V_gradB = 0;

%
U_uu = 10;      % U_up_up
U_dd = 10;      % U_dn_dn
U_ud = 10.0;      % U_up_dn

J_up = 1;
J_dn = 1;

staG0 = 0. * U;

mu_up_lt = -2*V_gradB*unit_tilT*(1:L) + mod(1:L,2)*staG0;
mu_dn_lt = V_gradB*unit_tilT*(1:L) + mod(1:L,2)*staG0;

BDC = 'obc';   % 'obc' - open boundary condition; 'pbs' - periodic case
% BDC = 'pbc';   % 'obc' - open boundary condition; 'pbs' - periodic case


%% ramp time list
nt = 101;
T = 30 * 1e-03;
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
fprintf('Elapsed time is %.6f seconds.\n',tEnd)

ham_all = hamiltonian_1d_bhm2C(basis,ham_elems,J_up,J_dn,...
    U_uu,U_dd,U_ud,mu_up_lt,mu_dn_lt);

ham = ham_all.ham;
fprintf('Number of nonzero hamiltonian elements: %d.\n\n',nnz(ham))


%% ground state info
k = 3;
[V,D] = eigs(ham,k,'sr');   % find k smallest real eigenvectors and eigenvalues

psi_gs = V(:,1);
energy_gs = D(1,1);
energy_gap = D(2,2) - D(1,1);

probl_gs = psi_gs .* conj(psi_gs);
figure('Color','w')
plot(probl_gs)
ax = gca;
ax.FontSize = 14;
xlim([1,ns])
xlabel('Basis index','FontSize',16)
ylabel('Probability','FontSize',16)


%%
density_up_Mt = [];
density_dn_Mt = [];

psic = phi_init;
probl = abs(conj(psic).*psic);

stat_nC = occupation_statistics_1d_bhm2C_Fcn(psic,basis);
density_up_Mt = cat(1,density_up_Mt,stat_nC.density_up_ltC);
density_dn_Mt = cat(1,density_dn_Mt,stat_nC.density_dn_ltC);

%
tStart = tic; 
for kk = 2:nt
    fprintf('Current process: %04d / %04d.\n',kk,nt)
    
    tSC0 = tic;
    % update time-dependent parameters   
    % ...
    % ...
    ham_cur = hamiltonian_1d_bhm2C(basis,ham_elems,J_up,J_dn,...
        U_uu,U_dd,U_ud,mu_up_lt,mu_dn_lt);
    hamC = ham_cur.ham;
    tEC0 = toc(tSC0);
    fprintf('time for generate hamiltonian: %.6f seconds.\n',tEC0)
    
    tSC = tic;
    psic = expv(-1i*dt,hamC,psic,1e-8,50);
    tEC = toc(tSC);
    fprintf('time for evolution: %.6f seconds.\n',tEC)
    
    stat_nC = occupation_statistics_1d_bhm2C_Fcn(psic,basis);
    density_up_Mt = cat(1,density_up_Mt,stat_nC.density_up_ltC);
    density_dn_Mt = cat(1,density_dn_Mt,stat_nC.density_dn_ltC);
    
    fprintf('\n')
end
tEnd = toc(tStart);
fprintf('Elapsed time is %.6f seconds.\n',tEnd)


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


%%
% prob_lt = psic .* conj(psic);
% sort(prob_lt)
