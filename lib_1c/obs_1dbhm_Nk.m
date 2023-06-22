function [s] = obs_1dbhm_Nk(psi,basis)
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
ql = linspace(-1*pi,1*pi,nl);


%% obs_adag_a
adag_a = struct();
note_mat = [];
vals_mat = [];

tic
for ia = 1:L
    for ib = 1:L
        op_c = operator_adag_a(basis,ia,ib);
        field_cur = ['site_',num2str(ia),'_',num2str(ib)];
        
        obs_val = psi'*(op_c*psi);
        adag_a.(field_cur) = obs_val;

        note_mat = cat(1,note_mat,[ia,ib]);
        vals_mat = cat(1,vals_mat,obs_val);

    end
end
toc


%%
nk_vals = zeros(1,nl);

tic
for kk = 1:nl
    qc = ql(kk);

    nk_vals(kk) = ...
        sum(exp(1i*qc*(note_mat(:,1)-note_mat(:,2))) ...
        .* vals_mat);
end
toc

% tic
% for kk = 1:nl
%     % if mod(kk,10) == 1
%     %     fprintf('Current: %d/%d\n',kk,nl)
%     % end
%     qc = ql(kk);
%     for ia = 1:L
%         for ib = 1:L
%             field_cur = ['site_',num2str(ia),'_',num2str(ib)];
%             nk_vals(kk) = nk_vals(kk) ...
%                 + exp(-1i*qc*(ia-ib))*adag_a.(field_cur);
% 
%         end
%     end
% end
% toc

nk_vals = abs(nk_vals);
Q = trapz(ql,nk_vals);
nk_vals = nk_vals/Q*basis.Nb;


s.ql = ql;
s.nk_vals = nk_vals;


% return
%%
figure('Color','w')
hold on
plot(ql,abs(nk_vals))
plot(ql,real(nk_vals))
plot(ql,imag(nk_vals))
hold off

