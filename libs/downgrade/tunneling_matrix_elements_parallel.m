function [elems_Jup,elems_Jdn] = tunneling_matrix_elements_parallel(basis,BDC,n_blocks)


if nargin < 3
    n_blocks = 800;
end
    

%%
elems_Jup = struct();
elems_Jdn = struct();


%%
n_bs = basis.n_bs;
L = basis.L;


%%
n_tasks = ceil(n_bs/n_blocks);
    
task_lop_st = zeros(n_tasks,1);
task_lop_ed = zeros(n_tasks,1);
task_lop_st(:) = (1:n_blocks:n_bs)';
task_lop_ed(1:n_tasks-1) = (n_blocks:n_blocks:n_bs-1)';
task_lop_ed(n_tasks) = n_bs;

s_taskAll = struct();
for kk = 1:n_tasks
    s_taskAll(kk).Jelems_infoA = [];
end

parfor kk = 1:n_tasks
    n_st = task_lop_st(kk);
    n_ed = task_lop_ed(kk);
    s_taskAll(kk).Jelems_infoA = ...
        tunneling_generate_subFcn(basis,BDC,n_st,n_ed);

end


%% index generate
idx_up_i_all = struct();
idx_up_j_all = struct();
idx_up_k_all = struct();
idx_dn_i_all = struct();
idx_dn_j_all = struct();
idx_dn_k_all = struct();
idx_up_ilt = [];
idx_up_jlt = [];
idx_up_klt = [];
idx_dn_ilt = [];
idx_dn_jlt = [];
idx_dn_klt = [];


switch BDC
    case 'obc'
        for kk = 1:L-1
            field_cur = ['site_',num2str(kk),'_',num2str(kk+1)];
            idx_up_i_all.(field_cur) = [];
            idx_up_j_all.(field_cur) = [];
            idx_up_k_all.(field_cur) = [];
            
            idx_dn_i_all.(field_cur) = [];
            idx_dn_j_all.(field_cur) = [];
            idx_dn_k_all.(field_cur) = [];

        end
        
    case 'pbc'
        for kk = 1:L
            field_cur = ['site_',num2str(kk),'_',num2str(mod(kk,L)+1)];
            idx_up_i_all.(field_cur) = [];
            idx_up_j_all.(field_cur) = [];
            idx_up_k_all.(field_cur) = [];
            
            idx_dn_i_all.(field_cur) = [];
            idx_dn_j_all.(field_cur) = [];
            idx_dn_k_all.(field_cur) = [];
            
        end
        
    otherwise
        error('Invalid boundary condition!')

end

fields = fieldnames(idx_up_i_all);
nfs = numel(fields);

%
for kk = 1:n_tasks
    sC = s_taskAll(kk).Jelems_infoA;
    
    for jj = 1:nfs
        fieldCur = fields{jj};
        
        idx_up_i_all.(fieldCur) = cat(1,idx_up_i_all.(fieldCur),...
            sC.idx_up_i_all.(fieldCur));
        idx_up_j_all.(fieldCur) = cat(1,idx_up_j_all.(fieldCur),...
            sC.idx_up_j_all.(fieldCur));
        idx_up_k_all.(fieldCur) = cat(1,idx_up_k_all.(fieldCur),...
            sC.idx_up_k_all.(fieldCur));
        
        idx_dn_i_all.(fieldCur) = cat(1,idx_dn_i_all.(fieldCur),...
            sC.idx_dn_i_all.(fieldCur));
        idx_dn_j_all.(fieldCur) = cat(1,idx_dn_j_all.(fieldCur),...
            sC.idx_dn_j_all.(fieldCur));
        idx_dn_k_all.(fieldCur) = cat(1,idx_dn_k_all.(fieldCur),...
            sC.idx_dn_k_all.(fieldCur));

    end
    
    idx_up_ilt = cat(1,idx_up_ilt,sC.idx_up_ilt);
    idx_up_jlt = cat(1,idx_up_jlt,sC.idx_up_jlt);
    idx_up_klt = cat(1,idx_up_klt,sC.idx_up_klt);
    idx_dn_ilt = cat(1,idx_dn_ilt,sC.idx_dn_ilt);
    idx_dn_jlt = cat(1,idx_dn_jlt,sC.idx_dn_jlt);
    idx_dn_klt = cat(1,idx_dn_klt,sC.idx_dn_klt);
    
end


%%
switch BDC
    case 'obc'
        for kk = 1:L-1
            field_cur = ['site_',num2str(kk),'_',num2str(kk+1)];
            
            i_up_cur = idx_up_i_all.(field_cur);
            j_up_cur = idx_up_j_all.(field_cur);
            v_up_cur = idx_up_k_all.(field_cur);
            
            i_dn_cur = idx_dn_i_all.(field_cur);
            j_dn_cur = idx_dn_j_all.(field_cur);
            v_dn_cur = idx_dn_k_all.(field_cur);
            
            elems_Jup.(field_cur) = ...
                sparse(i_up_cur,j_up_cur,v_up_cur,...
                n_bs,n_bs);
            elems_Jdn.(field_cur) = ...
                sparse(i_dn_cur,j_dn_cur,v_dn_cur,...
                n_bs,n_bs);
            
        end
        
    case 'pbc'
        for kk = 1:L
            field_cur = ['site_',num2str(kk),'_',num2str(mod(kk,L)+1)];
            
            i_up_cur = idx_up_i_all.(field_cur);
            j_up_cur = idx_up_j_all.(field_cur);
            v_up_cur = idx_up_k_all.(field_cur);
            
            i_dn_cur = idx_dn_i_all.(field_cur);
            j_dn_cur = idx_dn_j_all.(field_cur);
            v_dn_cur = idx_dn_k_all.(field_cur);
            
            elems_Jup.(field_cur) = ...
                sparse(i_up_cur,j_up_cur,v_up_cur,...
                n_bs,n_bs);
            elems_Jdn.(field_cur) = ...
                sparse(i_dn_cur,j_dn_cur,v_dn_cur,...
                n_bs,n_bs);
            
        end  
        
end
    
%
elems_Jup.sumAllSite = sparse(idx_up_ilt,idx_up_jlt,...
    idx_up_klt,n_bs,n_bs);
elems_Jdn.sumAllSite = sparse(idx_dn_ilt,idx_dn_jlt,...
    idx_dn_klt,n_bs,n_bs);

