% set relevance for intervals to nan that do not exist in a certain class
% but have been detected by the algorithm erroneously
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

function relv = setRelvToNan(relv, interv, labels)
% get metrics, calculated from relevance in intervals
fields_relv = string( fieldnames(relv) ); 

% reject 'RR' and 'filtRate'
fields_relv(ismember(fields_relv,{'RR','filtRate'})) = [];

% get index for P waves and F waves in intervals
mask_f  = ismember(interv,{'F'});
mask_p  = ismember(interv,{'P','PQ'});
idx_f   = find(mask_f);
idx_p   = find(mask_p);
mask_af = labels==1;
for f = fields_relv'
    % extract data from field
    temp = relv.(f);

    % set to nan
    if isa(temp, 'double')
        temp(mask_f,~mask_af,:) = 0;
        temp(mask_p,mask_af,:)  = 0;
    elseif isa(temp, 'cell')
        temp(mask_f,~mask_af)   = {[]};
        temp(mask_p,mask_af)    = {[]};
    end
  
    % give back to field
    relv.(f) = temp;
end
end