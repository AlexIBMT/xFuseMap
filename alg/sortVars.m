% sort the values of two related vectors
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

function idx = sortVars(in,var1,var2)
if size(var1,1) < size(var1,2)
    var1 = var1';
end
if size(var2,2) < size(var2,1)
    var2 = var2';
end

in2     = join(in,"_",2);
combs   = reshape(var1 + "_" + var2,[],1);
NC      = numel(combs);

idx     = zeros(NC,1);
ws      = 1;
for nc = 1:NC
    temp        = find(in2==combs(nc));
    we          = ws + numel(temp) - 1;
    idx(ws:we)  = temp;
    ws          = we + 1;
end
end