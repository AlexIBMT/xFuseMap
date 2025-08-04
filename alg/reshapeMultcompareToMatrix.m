% reshape results from multcompare to a matrix
% ----------------------------------------------------------------------- %
% For details see xFuseMap® main publication:
% A. Hammer et al., ‘Fusion of automatically learned rhythm and morphology
%   features matches diagnostic criteria and enhances AI explainability’, 
%   npj Artificial Intelligence, vol. 1, 2025,
%   doi: doi.org/10.1038/s44387-025-00022-w.
%
% Copyright
% Alexander Hammer
% Institute of Biomedical Engineering
% TU Dresden
% 01307 Dresden, Germany
%
% Version 1.0, Dresden 31.07.2025
% ----------------------------------------------------------------------- %

function [d,p,comb] = reshapeMultcompareToMatrix(in)
% INPUT ----------------------------------------------------------------- %
% in: multcompare matrix 'c' ('Group','Control Group','Lower Limit','Difference','Upper Limit','p-values')
% OUTPUT ---------------------------------------------------------------- %
% d: estimated difference between groups
% p: p values
% comb: combined [d1;p1;d2;p2;...dN;pN]
% ----------------------------------------------------------------------- %

% create empty arrays for d and p
g       = unique( in(:,1:2) );
ng      = numel(g);
[d,p]   = deal( nan(ng) );
comb    = nan(2*ng,ng);

% write values into d and p
nl      = size(in,1); % number of lines
for l = 1:nl
    i1 = in(l,1);
    i2 = in(l,2);
    d(i1,i2) = in(l,4);
    p(i1,i2) = in(l,end);
end

% combine d and p in comb
comb(1:2:2*ng,:) = d;
comb(2:2:2*ng,:) = p;