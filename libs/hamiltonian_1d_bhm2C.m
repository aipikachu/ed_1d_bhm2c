function [ham_mat] = hamiltonian_1d_bhm2C(basis,ham_elems,...
    J_up,J_dn,U_uu,U_dd,U_ud,mu_up_lt,mu_dn_lt)


%%
ham_mat = struct();


%%
n_bs = basis.n_bs;
L = basis.L;

BDC = ham_elems.BDC;


%% hamiltonian generation
ham = sparse(n_bs,n_bs);


% 01_1: Onsite interaction, (up,up);
len_U_uu = length(U_uu);
if len_U_uu == 1
    ham = ham + U_uu * ham_elems.ham_elems_U_uu.sumAllSite;
elseif len_U_uu == L
    for kk = 1:L
        field_cur = ['site_',num2str(kk)];
        ham = ham + U_uu(kk) * ham_elems.ham_elems_U_uu.(field_cur);
    end
else
    error('Error! Invalid input parameter "U_uu"!')
end


% 01_2: Onsite interaction, (dn,dn);
len_U_dd = length(U_dd);
if len_U_dd == 1
    ham = ham + U_dd * ham_elems.ham_elems_U_dd.sumAllSite;
elseif len_U_dd == L
    for kk = 1:L
        field_cur = ['site_',num2str(kk)];
        ham = ham + U_dd(kk) * ham_elems.ham_elems_U_dd.(field_cur);
    end
else
    error('Error! Invalid input parameter "U_dd"!')
end


% 01_3: Onsite interaction, (up,dn);
len_U_ud = length(U_ud);
if len_U_ud == 1
    ham = ham + U_ud * ham_elems.ham_elems_U_ud.sumAllSite;
elseif len_U_ud == L
    for kk = 1:L
        field_cur = ['site_',num2str(kk)];
        ham = ham + U_ud(kk) * ham_elems.ham_elems_U_ud.(field_cur);
    end
else
    error('Error! Invalid input parameter "U_ud"!')
end


% 02_1: chemical potential, |up>
len_mu_up = length(mu_up_lt);
if len_mu_up == 1
    ham = ham + mu_up_lt * ham_elems.ham_elems_mu_up.sumAllSite;
elseif len_mu_up == L
    for kk = 1:L
        field_cur = ['site_',num2str(kk)];
        ham = ham + mu_up_lt(kk) * ham_elems.ham_elems_mu_up.(field_cur);
    end
else
    error('Error! Invalid input parameter "mu_up_lt"!')
end


% 02_2: chemical potential, |dn>
len_mu_dn = length(mu_dn_lt);
if len_mu_dn == 1
    ham = ham + mu_dn_lt * ham_elems.ham_elems_mu_dn.sumAllSite;
elseif len_mu_dn == L
    for kk = 1:L
        field_cur = ['site_',num2str(kk)];
        ham = ham + mu_dn_lt(kk) * ham_elems.ham_elems_mu_dn.(field_cur);
    end
else
    error('Error! Invalid input parameter "len_mu_dn"!')
end


% 03: tunneling, |up> and |dn>
len_J_up = length(J_up);
len_J_dn = length(J_dn);

switch BDC
    case 'obc'
        if len_J_up == 1
            if ~sum(size(ham_elems.ham_elems_J_up.sumAllSite) == [0,0])
                ham = ham + J_up * ham_elems.ham_elems_J_up.sumAllSite; 
            end        
        elseif len_J_up == (L-1)
            for kk = 1:L-1
                field_cur = ['site_',num2str(kk),'_',num2str(kk+1)];
                if sum(size(ham_elems.ham_elems_J_up.(field_cur)) == [0,0])
                    continue
                end
                ham = ham + J_up(kk) * ham_elems.ham_elems_J_up.(field_cur);
            end
        else
            error('Error! Invalid input parameter "J_up"!')
        end
        
        if len_J_dn == 1
            if ~sum(size(ham_elems.ham_elems_J_dn.sumAllSite) == [0,0])
                ham = ham + J_dn * ham_elems.ham_elems_J_dn.sumAllSite;
            end
        elseif len_J_dn == (L-1)
            for kk = 1:L-1
                field_cur = ['site_',num2str(kk),'_',num2str(kk+1)];
                if sum(size(ham_elems.ham_elems_J_dn.(field_cur)) == [0,0])
                    continue
                end
                ham = ham + J_dn(kk) * ham_elems.ham_elems_J_dn.(field_cur);
            end
        else
            error('Error! Invalid input parameter "J_dn"!')
        end
        
    case 'pbc'
        if len_J_up == 1
            if ~sum(size(ham_elems.ham_elems_J_up.sumAllSite) == [0,0])
                ham = ham + J_up * ham_elems.ham_elems_J_up.sumAllSite; 
            end            
        elseif len_J_up == L
            for kk = 1:L
                field_cur = ['site_',num2str(kk),'_',num2str(mod(kk,L)+1)];
                if sum(size(ham_elems.ham_elems_J_up.(field_cur)) == [0,0])
                    continue
                end
                ham = ham + J_up(kk) * ham_elems.ham_elems_J_up.(field_cur);
            end
        else
            error('Error! Invalid input parameter "J_up"!')
        end
        
        if len_J_dn == 1
            if ~sum(size(ham_elems.ham_elems_J_dn.sumAllSite) == [0,0])
                ham = ham + J_dn * ham_elems.ham_elems_J_dn.sumAllSite;
            end
        elseif len_J_dn == L
            for kk = 1:L
                field_cur = ['site_',num2str(kk),'_',num2str(mod(kk,L)+1)];
                if sum(size(ham_elems.ham_elems_J_dn.(field_cur)) == [0,0])
                    continue
                end
                ham = ham + J_dn(kk) * ham_elems.ham_elems_J_dn.(field_cur);
            end
        else
            error('Error! Invalid input parameter "J_dn"!')
        end
        
        
end


%%
ham_mat.ham_elems = ham_elems;
ham_mat.ham = ham;


