function [ham_mat] = hamiltonian_1d_bhm(basis,ham_elems,J,U,mu_lt)


%%
ham_mat = struct();


%%
ns = basis.ns;
L = basis.L;

BDC = ham_elems.BDC;


%% hamiltonian generation
ham = sparse(ns,ns);


% 01: Onsite interaction;
len_U = length(U);
if len_U == 1
    ham = ham + U * ham_elems.ham_elems_U.sumAllSite;
elseif len_U == L
    for kk = 1:L
        field_cur = ['site_',num2str(kk)];
        ham = ham + U(kk) * ham_elems.ham_elems_U.(field_cur);
    end
else
    error('Error! Invalid input parameter "U"!')
end



% 02_1: chemical potential, |up>
len_mu = length(mu_lt);
if len_mu == 1
    ham = ham + mu_lt * ham_elems.ham_elems_mu.sumAllSite;
elseif len_mu == L
    for kk = 1:L
        field_cur = ['site_',num2str(kk)];
        ham = ham + mu_lt(kk) * ham_elems.ham_elems_mu.(field_cur);
    end
else
    error('Error! Invalid input parameter "mu_lt"!')
end



% 03: tunneling, 
len_J = length(J);

switch BDC
    case 'obc'
        if len_J == 1
            if ~sum(size(ham_elems.ham_elems_J.sumAllSite) == [0,0])
                ham = ham + J * ham_elems.ham_elems_J.sumAllSite; 
            end        
        elseif len_J == (L-1)
            for kk = 1:L-1
                field_cur = ['site_',num2str(kk),'_',num2str(kk+1)];
                if sum(size(ham_elems.ham_elems_J.(field_cur)) == [0,0])
                    continue
                end
                ham = ham + J(kk) * ham_elems.ham_elems_J.(field_cur);
            end
        else
            error('Error! Invalid input parameter "J"!')
        end
        
        
    case 'pbc'
        if len_J == 1
            if ~sum(size(ham_elems.ham_elems_J.sumAllSite) == [0,0])
                ham = ham + J * ham_elems.ham_elems_J.sumAllSite; 
            end            
        elseif len_J == L
            for kk = 1:L
                field_cur = ['site_',num2str(kk),'_',num2str(mod(kk,L)+1)];
                if sum(size(ham_elems.ham_elems_J.(field_cur)) == [0,0])
                    continue
                end
                ham = ham + J(kk) * ham_elems.ham_elems_J.(field_cur);
            end
        else
            error('Error! Invalid input parameter "J"!')
        end
        
end


%%
ham_mat.ham_elems = ham_elems;
ham_mat.ham = ham;


