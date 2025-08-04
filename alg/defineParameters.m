%%
% HEAD
% 
% 
% Code for publication:
% A. Hammer et al., ‘Fusion of automatically learned rhythm and morphology
%   features matches diagnostic criteria and enhances AI explainability’, 
%   npj Artificial Intelligence, vol. 1, 2025,
%   doi: doi.org/10.1038/s44387-025-00022-w.
%
% Copyright Alexander Hammer, TU Dresden, 31.07.2025

%% GENERAL PARAMETERS
sp.pxFuseMap     = 0; % 0: don't plot xFuseMaps; 1: plot xFuseMaps

sp.plotFalseClassifications = 0; % plot relevance for false classifications
                
sp.methods  = "dt";  % use only predefined methods ("dt")


%% F WAVE DETECTOR
fwdp.pp  = 41;              % min. peak prominence
fwdp.ps  = 0;               % min. peaks per second
fwdp.tp  = 2;               % threshold periodicity
fwdp.ti  = 0;               % threshold SD(IQR_f)/SD(IQR_t)
fwdp.m   = 'negative';      % mode
fwdp.exc = {'','QT'};       % segment exclusion


%% PLOT LABELS
% xFuseMap
lb.x    = 'Time / s';
lb.y    = 'Amplitude / a.u.';
lb.z    = 'rR / a.u.';
lb.sz   = 12;


%% COLOR FOR xFuseMap
clr.s       = 100;          % number of color splits  
clr.f0      = 0;            % factor of third color channel (combined relevance)
clr.f_lw    = 2.5;            % factor for line thickness in relation to relevance

clr.bs = [1 1 1];       % color center, standard: [0 0 0] = black
clr.m1 = [1 1 0];       % color model 1 > orange
clr.m2 = [0 1 1];       % color model 2 > blue
clr.bg = [1 1 1]*0.6;   % color background > gray
clr.f0 = 0.6;


%% RELEVANCE PER INTERVAL (Boxplots)
clr11       = [0.6350 0.0780 0.1840];
clr00       = [0 0.4470 0.7410];
clr20       = [0.4660 0.6740 0.1880];
clr10       = [0.8500 0.3250 0.0980];
clr01       = [0.3010 0.7450 0.9330];

res_box     = struct();
res_box.met = ["n","mean","std","median","perc25%","perc75%","iqr","min","max"];

ptype       = "Box";    % "Bar", "Box"
offset      = 0.1;