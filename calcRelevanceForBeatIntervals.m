% xFuseMap® main script
% plot xFuseMaps and do statistics on differences between the distribution
% of relevance values across intervals in the ECG on beat level
% ----------------------------------------------------------------------- %
% For details see main publication:
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


%% PREPARATIONS
% add algorithm directory to matlab path
addpath(genpath("alg"))

% load parameter set
defineParameters

% define paths
pth.dat     = fullfile(pwd,'data');         % data
pth.res     = fullfile(pwd,'results');      % results
if ~isfolder(pth.res)
    mkdir(pth.res)
end
pth.plt     = fullfile(pth.res,"plots");    % plots
if ~isfolder(pth.plt)
    mkdir(pth.plt)
end


%% LOAD DATA
load( fullfile(pth.dat,'xFuseMap_test_full.mat') );

% number of ECG recordings
nn  = numel(rec);

% classes
classes         = unique(labels, 'sorted');
classes_full    = unique(labels_full, 'sorted');


%% model names
models  = split(vn_pred,'_');
models  = unique(models(:,1));

% number of models
nm      = numel(models);


%% xFuseMap
% create label information for subtitle
st_lbl  = "Label: "+replace( string(labels_full),...
    {'0','1','2'},...
    {'n-AF (O)','AF','n-AF (NSR)'}) ;

% get decision certainty for majority class per model
[cert,dec]  = deal( zeros(nn,nm) );
for m = 1:nm
    mask.model      = contains(vn_pred,models(m));

    class           = extractAfter(vn_pred(mask.model),'_class_');
    
    [cert(:,m),dt]  = max(pred_weight(mask.model,:),[],1); % certainty
    dec(:,m)        = str2double( class(dt) );               % decision
end
clear('dt')

% iterate over recordings to create xFuseMap plots
if sp.pxFuseMap
    for n = 1:nn
        % load data (amplitude and relevance values [sig, relv])
        load(fullfile(pth.dat,'rec',rec(n)+".mat"));
    
        % time axis
        x = (1:numel(sig))/meta.Fs(n);
    
        % relevance values
        z = zeros(numel(sig),nm);
        for m = 1:nm
            mask.pred       = contains(vn_relv,models(m)) & ...
                contains(vn_relv,"_class_"+string(dec(n,m))) & ...
                contains(vn_relv,"_"+sp.methods);
            z(:,m)          = relv(:,:,mask.pred);
        end
    
        % title and subtitle
        models_t =  replace(models',{'longterm','shortterm'},{'long-term','short-term'});
        lb.t    = "Record: "+rec(n);
        lb.st   = join(...
            [models_t+": "+replace(string(dec(n,:)),{'0','1'},{'n-AF','AF'})+...
            " ("+string(round(cert(n,:)*100,1))+"%)",st_lbl(n)],...
            "   -   ");
    
        % create plot
        fig = figure;
        set(fig,'Position',[1 49 1920 955/2.5])
        fig = colorbarPlot(fig, x, sig, z, ...
            clr.s, clr.m1, clr.bs, clr.m2, clr.bg, clr.f0, clr.f_lw, ...
            lb.x, lb.y, lb.z, models_t', lb.sz, lb.t, lb.st, 1);
    
        % save
        fnm_temp = rec(n)+"_label-"+string(labels_full(n))+"_class-"+join(string(dec(n,:)),"-");
        exportgraphics(fig,fullfile(pth.plt,fnm_temp+".pdf"),...
            'ContentType','vector')
    
        % close figure
        close(fig)
        clear('fig','fnm_temp','models_t','x','sig','z','relv')
    end
end


%% do ANOVA
relv_int2           = setRelvToNan(relv_int, int_labels, labels); % eliminate P and PQ for AF and F for n-AF
[anov,anov_fields]  = doANOVA(rec, labels, labels_full, pred_weight, vn_pred, vn_relv, relv_int2, int_labels, models, [1 0 2], "AF");

% reshape and merge results
anov_rs         = struct();
anov_rs_n       = numel(classes_full)*numel(int_labels);
mask_all        = true(anov_rs_n);
mask_upper      = triu(mask_all,1);
mask_lower      = tril(mask_all,-1);
for m = models'
    [anov_rs.("d_mean_"+m),anov_rs.("p_mean_"+m),anov_rs.("comb_mean_"+m)] = reshapeMultcompareToMatrix(anov.(m+"_relv_mean_dt").ph_c);
    [anov_rs.("d_std_"+m),anov_rs.("p_std_"+m),anov_rs.("comb_std_"+m)] = reshapeMultcompareToMatrix(anov.(m+"_relv_std_dt").ph_c);
    
    j = 1;
    for k  = ["d","p"]
        anov_rs.(k+"_"+m) = nan(anov_rs_n);
        anov_rs.(k+"_"+m)(mask_upper) = anov_rs.(k+"_mean_"+m)(mask_upper);
        anov_rs.(k+"_"+m)(mask_lower) = anov_rs.(k+"_std_"+m)(mask_upper)';

        % writetable
        anov_rs.(k+"_"+m) = array2table(anov_rs.(k+"_"+m),"VariableNames",anov.ph_names,"RowNames",anov.ph_names);
        writetable(anov_rs.(k+"_"+m), ...
            fullfile(pth.res,"res_anova.xlsx"),...
            'Sheet',k+"_"+m,...
            'WriteVariableNames',true,'WriteRowNames',true)
    end
end


%% plot relevance
% ... per interval
% ... per class (speaks for AF vs speaks against AF)
% ... per classifier (shortterm vs longterm)

% iterate over metrics
for fi = ["relv_mean","relv_std"]    
    % get relevance values
    vals_all    = relv_int2.(fi);
    vals        = struct();

    % create figure
    fig         = figure;
    set(fig,'Position',[1 49,1920 955])

    % get method name
    met = sp.methods;

    for m = 1:nm
        mod = models{m};
        
        subplot(numel(sp.methods),numel(models),m)
        
        for c = classes
            mask.("relv"+string(c)) = contains(vn_relv,met) & contains(vn_relv,mod) & contains(vn_relv,string(c));
            mask.("pred"+string(c)) = contains(vn_pred,mod) & contains(vn_pred,string(c));

            vals.("v"+string(c))    = vals_all(:,:,mask.("relv"+string(c)))';
            % if c==1
            %     vals.("v"+string(c))(:,contains(interv,'P')) = 0;
            % end
        end

        % get relevance values for each class/label combination
        vals00  = vals.v0(labels_full==0 & pred(mask.pred0,:)>pred(mask.pred1,:),:);
        vals11  = vals.v1(labels_full==1 & pred(mask.pred1,:)>pred(mask.pred0,:),:);
        vals20  = vals.v0(labels_full==2 & pred(mask.pred0,:)>pred(mask.pred1,:),:);
        vals01  = vals.v1(labels_full==0 & pred(mask.pred1,:)>pred(mask.pred0,:),:);
        vals10  = vals.v0(labels_full==1 & pred(mask.pred0,:)>pred(mask.pred1,:),:);

        % reject outliers
        vals00(isoutlier(vals00,"median")) = nan;
        vals11(isoutlier(vals11,"median")) = nan;
        vals20(isoutlier(vals20,"median")) = nan;
        vals01(isoutlier(vals01,"median")) = nan;
        vals10(isoutlier(vals10,"median")) = nan;

        % calculate metrics
        res_box.(fi+"_"+met+"_"+mod+"_AF")    = calculateBoxMetrics(vals11,res_box.met);
        res_box.(fi+"_"+met+"_"+mod+"_O")     = calculateBoxMetrics(vals00,res_box.met);
        res_box.(fi+"_"+met+"_"+mod+"_NSR")   = calculateBoxMetrics(vals20,res_box.met);
        res_box.(fi+"_"+met+"_"+mod+"_ALL")   = calculateBoxMetrics([vals11;vals00;vals20],res_box.met);

        x = 1:numel(int_labels);
        % plot
        boxplot(vals11,'positions',x-0.35+offset,'width',0.2,'color',clr11,'PlotStyle','compact','Symbol','.','Jitter',0.15)
        hold on
        boxplot(vals00,'positions',x-0.1+offset,'width',0.2,'color',clr00,'PlotStyle','compact','Symbol','.','Jitter',0.15)
        boxplot(vals20,'positions',x+0.25,'width',0.2,'color',clr20,'PlotStyle','compact','Symbol','.','Jitter',0.15)
        xnew = reshape([-0.35+offset;-0.1+offset;0.25] + x,[],1)';
        ymin = reshape([min(vals11,[],1);min(vals00,[],1);min(vals20,[],1)],[],1);
        ymax = reshape([max(vals11,[],1);max(vals00,[],1);max(vals20,[],1)],[],1);
        highlightSigNiv2(fig,xnew,ymin,ymax,anov.(mod+"_"+fi+"_"+met).ph_cc,2,"",0.01,[])

        % prepare legend
        a = get(get(gca,'children'),'children');   % Get the handles of all the objects
        for idx_a = find(~cellfun(@isempty,a))'
            t = get(a{idx_a},'tag');          % List the names of all the objects 
            idx_elmnts=strcmpi(t,'box');      % Find Box objects
            elmnts=a{idx_a}(idx_elmnts);      % Get the children you need
            set(elmnts,'linewidth',8);
            idx_elmnts=strcmpi(t,'whisker');  % Find whisker objects
            elmnts=a{idx_a}(idx_elmnts);      % Get the children you need
            set(elmnts,'linewidth',1.25);
            idx_elmnts=strcmpi(t,'outliers'); % Find outlier objects
            elmnts=a{idx_a}(idx_elmnts);      % Get the children you need
            set(elmnts,'MarkerSize',5);
        end
        hold off
        
        % add labels
        xticks(x)
        xticklabels(string(int_labels))
        xlabel('Interval')
        ylabel(lb.z)
        if strcmp(fi,'relv_mean')
            ylim([-0.05 1.1])
        elseif strcmp(fi,'relv_std')
            ylim([-0.01 0.29])
        else
            ylim padded % tight, padded, tickaligned
        end
        xlim([0 numel(int_labels)]+0.5)
        title(replace(fi,"_"," ") + " - " + mod + " - " + met)

        % get figure elements
        elements    = findall(gca,'Tag',"Box");
        elements    = elements(sort(numel(sp.methods)*size(vals_all,1):size(vals_all,1):numel(elements),'descend'));
        
        % add legend
        lg = legend(elements,...
            ["AF    (n="+string(size(vals11,1))+")",...
            "O      (n="+string(size(vals00,1))+")",...
            "NSR (n="+string(size(vals20,1))+")"]);

        % set color and position
        set(lg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.7]));
        set(gca,'Position',[0.1+m/2-0.5 0.1 0.38 0.8],'FontSize',12)
    end

    % save figure
    exportgraphics(fig,fullfile(pth.res,"box_"+fi+".pdf"),...
       'ContentType','vector',...
       'Append',true)

    close(fig)
    clear fig
end

%% save
save(fullfile(pth.res,"res_all.mat"),"models","anov","res_box")