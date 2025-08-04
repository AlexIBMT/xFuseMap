% highlight significance levels
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

function highlightSigNiv2(fig,x,ymin,ymax,p_stat,withEndMark,ns,minSigLevel,rejectFields)
% INPUT ----------------------------------------------------------------- %
%   fig         figure
%   x           x values        (... x ...)     [double]
%   y           y values        (... x ...)     [double]
%   p_stat      p statistics    (... x ...)     [double]
%   withEndMark 0: significance bar without end marks
%               1: ... with end marks
%               2: ... with lines from significance level to maximum value
%   ns          "":   show non-significant relations as well as significant relations
%               "only":   only highlight non-significant relations
%               "ignore": ignore non-significant relations
%   minSigLevel minimum significance level (0...1) [double]
% OUTPUT ---------------------------------------------------------------- %
%   -
% ----------------------------------------------------------------------- %
y = 1;
cl  = unique(x,'sorted');
if size(cl,1) < size(cl,2)
    cl = cl';
end

if ~exist('minSigLevel','var') || isempty(minSigLevel)
    minSigLevel = 0.05;
end

% sort array (short connections > long connections)
idxDiff     = abs(diff(p_stat(:,1:2),1,2));
[~,idxNew]  = sort(idxDiff,1,'ascend');
p_stat      = p_stat(idxNew,:);

% create array with significancy values
Na      = size(p_stat,1);
sig_niv = cell(Na,1);
if ~strcmp(ns,'only')
    if minSigLevel==0.05
        sig_niv(p_stat(:,end)<=0.05)        = {"*"};
        sig_niv(p_stat(:,end)<0.01)         = {"**"};
        sig_niv(p_stat(:,end)<0.001)        = {"***"};
    else
        sig_niv(p_stat(:,end)<=minSigLevel) = {"*"};
    end
end
if ~strcmp(ns,'ignore')
    sig_niv(p_stat(:,end)>minSigLevel)      = {"ns"};
end

% plot
hold on
levels          = nan(Na,1);
sig_nivEmpty    = ~cellfun(@isempty,sig_niv);
for a = 1:Na
    if ~isempty(sig_niv{a}) && ~any(ismember(p_stat(a,1:2),rejectFields))
        % get number of significance lines, that cross the current line
        clo = cl(p_stat(1:a,1:2));
        if size(clo,1) ~= a
            clo = clo';
        end
%         k   = max(levels(~cellfun(@isempty,sig_niv(1:a)) & ...
        if sum(sig_nivEmpty(1:a))==1
            k = 1;
        else
            m = 1:a-1;
            levels_used = unique(levels(sig_nivEmpty(m) & ...
                ((clo(a,1) > clo(m,1) & clo(a,1) < clo(m,2)) | ...          % left coordinate is inside another couple comparison
                (clo(a,2) > clo(m,1) & clo(a,2) < clo(m,2)) | ...           % left coordinate is inside another couple comparison
                (clo(a,1) >= clo(m,1) & clo(a,2) <= clo(m,2)) | ...         % couple comparison is completely inside another couple comparison
                (clo(m,1) >= clo(a,1) & clo(m,2) <= clo(a,2)))));           % another couple comparison is complete inside this couple comparison
            levels_all = 1:max(levels)+1;
            k = min(levels_all(~ismember(levels_all,levels_used)));
        end
        levels(a) = k;

        % plot
        sig_x   = linspace(cl(p_stat(a,1)),cl(p_stat(a,2)),2);
        sig_yp  = ones(1,numel(sig_x)).*(max(ymax)+diff([min(ymin),max(ymax)])*0.05*k);
        if strcmp(sig_niv{a},"ns")
            fyt = 0.3;
        else
            fyt = 0.15;
        end
        sig_yt  = max(ymax)+diff([min(ymin),max(ymax)])*0.05*(k+fyt);
        plot(sig_x, sig_yp, 'k-')
        if withEndMark >= 1 % line boundary
            fEM = diff([min(ymin),max(ymax)])*0.008;
%             xl = xline(sig_x,'Color',[1 1 1]*0.7);
%             uistack(xl,'bottom')
            plot(repmat([sig_x(1) sig_x(end)],2,1),sig_yp(1:2)-[0;fEM],'k-');
        end
        if withEndMark == 2 % 
            plot(repmat([sig_x(1) sig_x(end)],2,1),[ymax(p_stat(a,[1 2]))';sig_yp(1:2)],'k:')
        end
        text(mean(sig_x), sig_yt, sig_niv{a}, 'HorizontalAlignment', 'center')
        ylim padded
    end
end
hold off