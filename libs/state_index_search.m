function [idx_st] = state_index_search(st_up,st_dn,basis)


%%
nMax = basis.nMax;
idxlt = basis.idxstatel;


%%
idxIntl_up = base2dec(strrep(num2str(st_up),' ',''),nMax+1);
idxIntl_dn = base2dec(strrep(num2str(st_dn),' ',''),nMax+1);


%%
[~,idx_st] = ismember([idxIntl_up,idxIntl_dn],idxlt,'rows');

