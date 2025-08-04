% script for calculating the analysis of variance for different metrics
% that summarize the distribution of relevance values across ECG recordings
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

function [stat,fs] = doANOVA(names,labels_used,labels_full,pred,vn_pred,vn_relv,relv_int,interv,models,classes,patho)
vn_pred = strrep(vn_pred,"-","");
vn_relv = strrep(vn_relv,"-","");
models  = strrep(models,"-","");

if size(labels_used,2) > size(labels_used,1)
    labels_used = labels_used';
end
if size(labels_full,2) > size(labels_full,1)
    labels_full = labels_full';
end

% models
xAIs    = split(string(vn_relv),'_');
xAIs    = unique(xAIs(:,end));
nX      = numel(xAIs);

% classes
if ~exist('classes','var') || isempty(classes)
    classes = sort(unique(labels_full),'ascend');
end
nC      = numel(classes);

% 
nN      = numel(names);

%
hotEnc  = all( contains(vn_pred,"class_") );

% get prediction per model & mask for correct classification
nMo     = numel(models);
if hotEnc
    ppMod   = zeros(numel(labels_full),nMo);
    for mo = 1:nMo
        mask.mo     = contains(vn_pred,models{mo});
        mask.mo0    = mask.mo & contains(vn_pred,"class_0");
        mask.mo1    = mask.mo & contains(vn_pred,"class_1");
        ppMod(pred(mask.mo1,:)>pred(mask.mo0,:),mo) = 1;
    end
    mask    = rmfield(mask,{'mo','mo0','mo1'});
else
    ppMod = pred;
end

% mask.fc (false classified)
mask.fc = labels_used ~= ppMod;

% metrics (mean, meanMean_wv, std, stdMean_wv) -> nMe = 4
metr    = {'relv_mean','relv_std','drRmean'};
metr    = metr(ismember(metr,fieldnames(relv_int)));
nMe     = numel(metr);

% intervals (p, pq, q, r, s, st, t, qt, tq, f) -> nI = 10
nI      = numel(interv);


%% create one array per model & metric
arr     = struct();
for x = 1:nX
    for mo = 1:nMo
        for me = 1:nMe
            % create field name
            nm          = models{mo}+"_"+metr{me}+"_"+xAIs(x);
    
            % mask for relevance data
            if hotEnc
                mask.relv   = repmat(vn_relv',nN,1)==models{mo}+"_class_"+string(ppMod(:,mo))+"_"+xAIs(x);
            else
                mask.relv   = repmat(vn_relv',nN,1)==models{mo}+"_"+xAIs(x);
            end
                
            % extract relevance
            temp        = zeros(nI,nN);
            for r = 1:numel(vn_relv)
                temp(:,mask.relv(:,r)) = relv_int.(metr{me})(:,mask.relv(:,r),r);
            end
    
            % remove false classified
            temp(:,mask.fc(:,mo)) = nan;

            % save with outliers
            arr.(nm+"_wOutl") = temp;

            % remove outliers
            for c = 1:nC
                mask.c              = find(labels_full == classes(c));
                mask.outl           = isoutlier(temp(:,mask.c),"median",2);
                mask.comb           = false(size(temp));
                mask.comb(:,mask.c) = mask.outl;

                temp(mask.comb)     = nan;
            end
    
            % save without outliers
            arr.(nm)    = temp;
        end
    end
end
clear('temp')


%% ANOVA
fs      = string( fields(arr) );    % fields
stat    = struct();                 % statistics

% fit labels & intervals
labels2 = string(reshape( repmat(labels_full',nI,1), [],1));
interv2 = reshape( repmat(string(interv'),1,nN), [],1);

% sort variables
idx     = sortVars([labels2 interv2],string(classes),string(interv));
labels2 = labels2(idx);
interv2 = interv2(idx);
for f = fs'
    % reshape
    dat = reshape( arr.(f), [],1);
    dat = dat(idx);
    
    % ANOVA
    [~,stat.(f).anova,stat.(f).stats] = anovan(...
        dat, {labels2,interv2},...
        'model',2,...
        'varnames',{'labels','intervals'},...
        'display','off');

    % check for normal distribution of residuals (p > 0.05 -> normal distribution)
    [stat.(f).ks_h, stat.(f).ks_p] = kstest(stat.(f).stats.resid);
    
    % post-hoc test
    [stat.(f).ph_c, stat.(f).ph_m, ~, stat.ph_names] = multcompare(...
        stat.(f).stats,...
        "Dimension",[1 2],...
        'display','off');
    
    % create matrix
    temp1   = split(stat.ph_names,[",","="]);
    temp1   = string(temp1(:,[2 4]));
    stat.(f) = createResTable(nI,nC,0,stat.(f),"ph",temp1);     % segment and class comparison
    stat.(f) = createResTable(nI,nC,1,stat.(f),"ph_sc",temp1);  % segment comparison
    stat.(f) = createResTable(nI,nC,2,stat.(f),"ph_cc",temp1);  % class comparison

    % reject fields
    if any(strcmp(patho,{'AF','AFlut'}))
        rjf     = find((ismember(temp1(:,1),["0";"2"]) & ismember(temp1(:,2),"F")) | ...
            (ismember(temp1(:,1),"1") & ismember(temp1(:,2),["P";"PQ"])));
    else
        rjf     = [];
    end
    rj_idx  = find(any(ismember(stat.(f).ph(:,[1 2]),rjf),2));
    stat.(f).ph(rj_idx,:)       = [];
    rj_idx  = find(any(ismember(stat.(f).ph_sc(:,[1 2]),rjf),2));
    stat.(f).ph_sc(rj_idx,:)    = [];
    rj_idx  = find(any(ismember(stat.(f).ph_cc(:,[1 2]),rjf),2));
    stat.(f).ph_cc(rj_idx,:)    = [];
end
end

function stat = createResTable(nI,nC,idx,stat,nm,temp1)
    idx_mat = zeros(nI*nC);
    for i = 1:nI*nC
        if idx > 0
            idx_mat(i,:)    = temp1(:,idx)==temp1(i,idx);
        else
            idx_mat(i,:)    = any(temp1==temp1(i,:),2);
        end
    end
    idx_mat = triu(idx_mat,1);
    
    [t_row,t_col] = find(idx_mat);
    tt = [t_row,t_col];
    tt = sortrows(tt,[1 2]);

    % new res table
    temp2 = ismember(join(string(stat.ph_c(:,[1 2])),',',2), join(string(tt),',',2));
    stat.(nm) = stat.ph_c(temp2,:);
end