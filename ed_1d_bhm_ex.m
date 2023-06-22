close all
clc

addpath(genpath('libs'))
addpath(genpath('lib_1c'))
addpath(genpath('src'))


%% init state
L = 10;
% st_init = [1 1 1 1 1 1];
st_init = ones(1,L);

if (length(st_init) ~= L)
    error('Error! Invalid input parameters!')
end

%
N_atom = sum(st_init);

% the maximum occupation atom number in each lattice site
nMax = 3; 


%% parameters setting
% Hubbard parameters
J = 200 * 2 *  pi;
U = 400 * 2 *  pi;

% gradient field
V_gradB = 0;

% chemical potential
mu0 = 0 * 2 * pi;
staG0 = 0.0 * U;
mu_lt = mu0 + V_gradB + mod(1:L,2)*staG0;

% boundary condition
% 'obc' - open boundary condition;
% 'pbc'  periodic case
BDC = 'obc';   
% BDC = 'pbc';   


%% basis generate
basis = boson_basis_1d(L,N_atom,nMax);
ns = basis.ns;
fprintf('Total basis number is %d.\n',ns)


%% state index search
idx_init_st = state_in_index(st_init,basis);

% basis.state(idx_init_st,:)
phi_init = zeros(ns,1);
phi_init(idx_init_st) = 1;


%% hamiltonian elements generate
tStart = tic; 
ham_elems = hamiltonian_1d_bhm_elements(basis,BDC);
tEnd = toc(tStart);
fprintf('Elapsed time is %.6f seconds.\n',tEnd)


%%
ham_cur = hamiltonian_1d_bhm(basis,ham_elems,J,U,mu_lt);

ham = ham_cur.ham;
fprintf('Number of nonzero hamiltonian elements: %d.\n',nnz(ham))


%% ground state info
k = 3;
% find k smallest real eigenvectors and eigenvalues
[V,D] = eigs(ham,k,'sr');   

psi_gs = V(:,1);
energy_gs = D(1,1);

probl_gs = psi_gs .* conj(psi_gs);
figure('Color','w')
plot(probl_gs)

% energy_gs/(probl_gs'*(full(ham)*probl_gs))

%% measurement
s_n = obs_1dbhm_N(psi_gs,basis);
s_nn = obs_1dbhm_NN(psi_gs,basis);


s_nk = obs_1dbhm_Nk(psi_gs,basis);
s_nknk = obs_1dbhm_NkNk(psi_gs,basis);

