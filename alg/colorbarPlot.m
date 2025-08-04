% plot amplitude values y over coordinates x in colors defined by z
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

function fig = colorbarPlot(fig,x,y,z,s,clr_1,clr_0,clr_2,clr_b,f0,f_lw,x_label,y_label,z_label,edges,fontsz,tit,subtit,plotMap)

if isempty(clr_0)
    clr_0 = [0 0 0];
end
clr_all = 1 - clr_0; 

%% transform input data
if size(x,1) > size(x,2)
    x = x';
end

if size(y,1) > size(y,2)
    y = y';
end

if size(z,1) > size(z,2)
    z = z';
end
J = size(z,1);

if size(x,2)>size(y,2) 
    x = x(1:size(y,2));
end

if J > 1
    z(z<0)  = 0;
end
x_3 = mkSegs(x);            % choose data points between values (x and y)
y_3 = mkSegs(y);
if max(abs(z),[],'all')>1
    z_n = z./max(abs(z),[],2);  % normalize to maximum
else
    z_n = z;
end

% replace nans by zeros
z_n(isnan(z_n)) = 0;


%% create matrix with colors
clrs = ones(s,s,3,4);
   
l    = zeros(3,1);
l(1) = find(clr_1 ~= mode(clr_1));
l(2) = find(clr_2 ~= mode(clr_2));
l(3) = find(~ismember([1 2 3], [l(1) l(2)]));

%
idx              = linspace(clr_0(1),clr_all(1),s);
clrs(:,:,l(1),1) = repmat(idx',1,s,numel(l(1)));
clrs(:,:,l(2),2) = repmat(idx,s,1,numel(l(2)));
clrs(:,:,l(3),3) = ones(s).*clr_0(l(3));
temp = ones(s)*f0.*idx;
if clr_all(1) < clr_0(1)
    clrs(:,:,l(3),4) = max(temp,flip(temp)');
else
    clrs(:,:,l(3),4) = min(temp,flip(temp)');
end

if f0 == 0
    clrs = min(clrs(:,:,:,1:3),[],4);
elseif f0 > 0 && mean(clr_1) < 0.5
    clrs = min(clrs(:,:,:,[1 2 4]),[],4);
else
    clrs = rot90(fliplr(min(clrs(:,:,:,[1 2 4]),[],4)),-1);
end


%% assign color to each sample
I = numel(y);
idx_clrs = zeros(I,J);
for j = 1:J
    [~,idx_clrs(:,j)]   = min(abs(idx-z_n(j,:)'),[],2);
end
if mean(clr_1)>0.5
    idx_clrs = abs(idx_clrs-(s+1));
end
clrs_fit = zeros(I,3);
for i = 1:I
    clrs_fit(i,:) = squeeze(clrs(idx_clrs(i,1),idx_clrs(i,2),:));
end
clrs_fit = clrs_fit + abs(min(min(clrs_fit,[],'all'),0));


%% plot
if numel(size(clrs))<=2 
    c = colorbar;
    colormap(clrs);
    set(gca,'YDir','normal')
else
    pos_old = get(gca,'Position');
    pos_fig = get(fig,'Position');
    pos_new = [pos_old(1)-0.05, pos_old(2)+0.05, pos_old(3), min(0.68,pos_old(4))];
    p4      = pos_old(3)*0.5;
    pos_map = [pos_new(1)+pos_new(3)*1.05, pos_old(2)+(pos_old(4)-p4)/2, p4*pos_fig(4)/pos_fig(3), p4];

    if plotMap
        subplot('Position',pos_map)
        image(clrs)
        set(gca,'FontSize',fontsz,'YDir','normal')
        xticks([1 size(clrs,2)])
        xticklabels([0 1])
        yticks([1 size(clrs,2)])
        yticklabels([0 1])
        ylabel(edges{1})
        xlabel(edges{2})
        subtitle(z_label)
    end
    
    subplot('Position',pos_new)
end
if ~isempty(edges) && numel(size(clrs))<=2
    c.Ticks = [-1 1];
    c.TickLabels = {edges{2},edges{1}};
    caxis([min(idx),max(idx)])
end
% hold on
% I = numel(x);
% for i = 1:I
%     if mean(clr_1)<0.5 || (mean(clr_1)>0.5 && f0>0)
%         w = 4*max(clrs_fit(i,:));
%     else
%         w = 4*(1-min(clrs_fit(i,:)));
%     end
%     plot(x_3(:,i),y_3(:,i),'Color',clrs_fit(i,:),'LineWidth',w+0.5);
% 
% 
% end
xvec    = [x(1:I-1);x(2:I)];
yvec    = [y(1:I-1);y(2:I)];
H       = plot(xvec, yvec);
lw      = max(abs(z_n),[],1)*f_lw + 0.7;
for i = 1:I-1
    H(i).LineWidth = lw(i);
end
colororder(clrs_fit(:,:))
% hold off

xlim tight
xlabel(x_label);
ylabel(y_label);
if numel(size(clrs))<=2
    c.Label.String = z_label;
    set(c,'FontSize',fontsz)
end

title(tit)
subtitle(subtit)

set(gca,'FontSize',fontsz)
if ~isempty(clr_b)
    set(gca,'color',clr_b)
end
grid on
end

%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = mkSegs(in)
    out = [nan in(2:end)-diff(in)/2;  in;  in(1:end-1)+diff(in)/2 nan];
end
% ----------------------------------------------------------------------- %

function out = linspace2(pl,pc,pr,s)
    sl = size(pl);
    sr = size(pr);
    pc = ones(max(sl,sr)).*pc;

    nc = numel(pc);
    if ~isempty(pl) && ~isempty(pr)
        out = zeros(2*s+1,nc);
    else
        out = zeros(s+1,nc);
    end

    [outl, outr] = deal([]);
    for i = 1:nc
        if ~isempty(pl)
            if pl(i) == pc(i)
                outl = ones(1,s+1)*pl(i);
            else
                outl = pl(i):(pc(i)-pl(i))/s:pc(i);
            end
        end
        if ~isempty(pr)
            if pr(i) == pc(i)
                outr = ones(1,s+1)*pr(i);
            else
                outr = pc(i):(pr(i)-pc(i))/s:pr(i);
            end
        end
        if ~isempty(pl) && ~isempty(pr)
            outr = outr(2:end);
        end
        out(:,i) = [outl outr];
    end
end

