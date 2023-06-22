function [s] = obs_1dbhm_NkNk(psi,basis)
%% [] = obs_1dbhm_Nk(psi_in,basis)
% state representation
% |psi> = \sum_i \alphi_i |\alphi_i>

%%
s = struct();


%%
L = basis.L;
% ns = basis.ns;
% nMax = basis.nMax;
% state = basis.state;


%%
nl = 101;
ql = linspace(-2*pi,2*pi,nl);


%% obs_adag_a
adag_a = struct();

tic
for ia = 1:L
    for ib = 1:L
        op_c = operator_adag_a(basis,ia,ib);
        field_cur = ['site_',num2str(ia),'_',num2str(ib)];

        adag_a.(field_cur) = psi'*(op_c*psi);

    end
end
toc

%% obs_adag_a_adag_a
ad_a_ad_a = struct();
note_mat = [];
vals_mat = [];

tic
for ia = 1:L
    for ib = 1:L
        op_ab = operator_adag_a(basis,ia,ib);
        for ic = 1:L
            for id = 1:L
                op_cd = operator_adag_a(basis,ic,id);
                field_cur = ['site_',num2str(ia),'_',...
                    num2str(ib),'_',num2str(ic),...
                    '_',num2str(id)];

                obs_val = psi'*(op_ab*(op_cd*psi));
                ad_a_ad_a.(field_cur) = obs_val;
                    
                note_mat = cat(1,note_mat,[ia,ib,ic,id]);
                vals_mat = cat(1,vals_mat,obs_val);
            end
        end
    end
end
toc

ad_a_ad_a.note_mat = note_mat;
ad_a_ad_a.vals_mat = vals_mat;


%%
nknk_vals = zeros(nl,nl);

tic
for iy = 1:nl
    for ix = 1:nl
        qx = ql(ix);
        qy = ql(iy);
        nknk_vals(iy,ix) = ...
            sum(exp(1i*qx*(note_mat(:,1)-note_mat(:,2))) ...
            .* exp(1i*qy*(note_mat(:,3)-note_mat(:,4))) ...
            .* vals_mat);

    end
end
toc


%%
figure('Color','w')
imagesc(real(nknk_vals))


s.ql = ql;
s.nknk_vals = nknk_vals;
s.ad_a_ad_a = ad_a_ad_a;

return

%%
