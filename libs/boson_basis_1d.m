function [basis] = boson_basis_1d(L,Nb,nMax)
%BOSON_
% boson_basis_1d(L,Nb,nMax), a draft of 
%  generate ...
%
%  INPUT: 
%        - 
%
%  OUTPUT:
%     basis -  A 
%
%  ---------------------------------
%   Author: Wei-Yong Zhang
%   Email: zhangwya@mail.ustc.edu.cn
%          zhangwya@gmail.com
%   Revision: 1.0
%   Create Date: Nov. 15, 2020
%   Last Update: Sep. 24, 2022
%  ---------------------------------


%%
basis = struct();


% All possible placements of internal dividers.                                       
dividers = nchoosek(1:(Nb+L-1),L-1);                                                  
ndividers = size(dividers,1);                                                        
% Add dividers at the beginning and end.                                              
b = cat(2,zeros(ndividers,1),dividers,(Nb+L)*ones(ndividers,1));                  
% Find distances between dividers.                                                    
c = diff(b,1,2)-1;
d = flip(c,1);

% find index 
e = sum(d <= ones(1,L)*nMax,2) == L;

state = d(e,:);
ns = size(state,1);

% idxIntl = zeros(ns,1);
% for kk = L:-1:1
%     idxIntl = idxIntl + state(:,kk) * (nMax+1)^(L-kk);
% end
idxIntl = state * ((nMax+1).^(L-1:-1:0))';


%%
basis.state = state;
basis.ns = ns;
basis.idxIntl = idxIntl;
basis.nMax = nMax;
basis.L = L;
basis.Nb = Nb;

