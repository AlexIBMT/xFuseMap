% script for calculating metrics to describe the results shown in boxplots
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

function out = calculateBoxMetrics(vals,vars)
    NV  = numel(vars);
    out = nan(NV,size(vals,2));
    for nv = 1:NV
        if strcmp(vars(nv),"iqr")
            out(nv,:) = iqr(vals,1);
        elseif strcmp(vars(nv),"max")
            out(nv,:) = max(vals,[],1);
        elseif strcmp(vars(nv),"mean")
            out(nv,:) = mean(vals,1,'omitnan');
        elseif strcmp(vars(nv),"median")
            out(nv,:) = median(vals,1,'omitnan');
        elseif strcmp(vars(nv),"min")
            out(nv,:) = min(vals,[],1);
        elseif strcmp(vars(nv),"n")
            out(nv,:) = sum(~isnan(vals));
        elseif strcmp(vars(nv),"perc25%")
            out(nv,:) = prctile(vals,25,1);
        elseif strcmp(vars(nv),"perc75%")
            out(nv,:) = prctile(vals,75,1);
        elseif strcmp(vars(nv),"std")
            out(nv,:) = std(vals,1,'omitnan');
        end
    end
end