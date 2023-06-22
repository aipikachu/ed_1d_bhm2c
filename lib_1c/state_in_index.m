function [idx_st] = state_in_index(st_in,basis)


%%
nMax = basis.nMax;
idxlt = basis.idxIntl;


%%
idxIntl_up = base2dec(strrep(num2str(st_in),' ',''),nMax+1);


%%
[~,idx_st] = ismember(idxIntl_up,idxlt);
