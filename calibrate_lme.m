% Mixed effects models of tree-ring stress index
n_min = 75;
syear = 1971;

% Number of vars for each environmental stress type
nTopo = 4;
nSoil = 7; % Subject to change
nTemp = 8;
nWater = 8;

% Names of variables
TopoLabs = {'ELEV','SLOPE','UAA','TWI'};
SoilLabs = {'CLAY','OC','PHCA','SAND','SILT','TN','VMC2'};
TempLabs = {'TMIN_SON','TMIN_DJF','TMIN_MAM','TMIN_JJA',...
	'TMAX_SON','TMAX_DJF','TMAX_MAM','TMAX_JJA'};
WaterLabs = {'WB_SON','WB_DJF','WB_MAM','WB_JJA',...
	'VPD_SON','VPD_DJF','VPD_MAM','VPD_JJA'};
mTempLabs = {'mTMIN_SON','mTMIN_DJF','mTMIN_MAM','mTMIN_JJA',...
	'mTMAX_SON','mTMAX_DJF','mTMAX_MAM','mTMAX_JJA'};
mWaterLabs = {'mWB_SON','mWB_DJF','mWB_MAM','mWB_JJA',...
	'mVPD_SON','mVPD_DJF','mVPD_MAM','mVPD_JJA'};
AllLabs = ['Stress','Year',TopoLabs SoilLabs TempLabs WaterLabs mTempLabs mWaterLabs,'Site','Species'];

p = length(AllLabs);

% month of previous year to start with
smonth = 9; % September

%% Load tree-ring stress estimates and put in table
load ./data/TREESI_wEnvFactors;

n = sum([TREESI.END] - syear + 1);

T = table('size',[n 2+nTopo+nSoil+nTemp*2+nWater*2+2], 'VariableTypes',...
    {'double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double',...
    'double','double','double','string','string'},...
    'VariableNames',AllLabs);
idx0 = 1;

for i = 1:length(TREESI)
    
    yr = TREESI(i).YEAR;
    ind = yr >= syear;
    idx = idx0:(idx0+sum(ind)-1);
    idx0 = idx0+sum(ind);
    
    T.Stress(idx) = TREESI(i).TREESI(ind);
    T.Year(idx) = yr(ind);
    T.ELEV(idx) = TREESI(i).ELEV;
    T.SLOPE(idx) = TREESI(i).SLOPE;
    T.UAA(idx) = TREESI(i).UAA;
    T.TWI(idx) = TREESI(i).TWI;
    T.CLAY(idx) = TREESI(i).GSDE.CLAY;
    T.OC(idx) = TREESI(i).GSDE.OC;
    T.PHCA(idx) = TREESI(i).GSDE.PHCA;
    T.SAND(idx) = TREESI(i).GSDE.SAND;
    T.SILT(idx) = TREESI(i).GSDE.SILT;
    T.TN(idx) = TREESI(i).GSDE.TN;
    T.VMC2(idx) = TREESI(i).GSDE.VMC2; % volumetric water content at -33 kPa

    T.Site(idx) = TREESI(i).SITE;
    T.Species(idx) = TREESI(i).SPECIES;
    
    [~, ind] = intersect(TREESI(i).PRISM.YEAR, yr(ind));
    
    % 3-month P-PET sums
    P = reshape((TREESI(i).PRISM.PPT - TREESI(i).PRISM.PET)', 1, []);
    windowSize = 3;
    b = ones(1,windowSize);
    a = 1;
    P = filter(b, a, P);
    P = reshape(P, 12, [])';
    T.WB_SON(idx) = P(ind-1, 11)- mean(P(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 11));
    T.WB_DJF(idx) = P(ind, 2)- mean(P(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 2));
    T.WB_MAM(idx) = P(ind, 5)- mean(P(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 5));
    T.WB_JJA(idx) = P(ind, 8)- mean(P(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 8));
    T.mWB_SON(idx) = mean(P(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 11));
    T.mWB_DJF(idx) = mean(P(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 2));
    T.mWB_MAM(idx) = mean(P(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 5));
    T.mWB_JJA(idx) = mean(P(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 8));
    
    % 3-month Tmin, Tmax, and VPD means
    Tmin = reshape(TREESI(i).PRISM.TMIN', 1, []);
    Tmax = reshape(TREESI(i).PRISM.TMAX', 1, []);
    VPD = reshape(TREESI(i).PRISM.VPD', 1, []);
    windowSize = 3;
    b = ones(1,windowSize)/windowSize;
    a = 1;
    Tmin = reshape(filter(b, a, Tmin), 12, [])';
    Tmax = reshape(filter(b, a, Tmax), 12, [])';
    VPD = reshape(filter(b, a, VPD), 12, [])';
    T.TMIN_SON(idx) = Tmin(ind-1, 11)- mean(Tmin(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 11));
    T.TMIN_DJF(idx) = Tmin(ind, 2)- mean(Tmin(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 2));
    T.TMIN_MAM(idx) = Tmin(ind, 5)- mean(Tmin(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 5));
    T.TMIN_JJA(idx) = Tmin(ind, 8)- mean(Tmin(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 8));
    T.mTMIN_SON(idx) = mean(Tmin(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 11));
    T.mTMIN_DJF(idx) = mean(Tmin(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 2));
    T.mTMIN_MAM(idx) = mean(Tmin(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 5));
    T.mTMIN_JJA(idx) = mean(Tmin(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 8));
    T.TMAX_SON(idx) = Tmax(ind-1, 11)- mean(Tmax(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 11));
    T.TMAX_DJF(idx) = Tmax(ind, 2)- mean(Tmax(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 2));
    T.TMAX_MAM(idx) = Tmax(ind, 5)- mean(Tmax(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 5));
    T.TMAX_JJA(idx) = Tmax(ind, 8)- mean(Tmax(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 8));
    T.mTMAX_SON(idx) = mean(Tmax(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 11));
    T.mTMAX_DJF(idx) = mean(Tmax(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 2));
    T.mTMAX_MAM(idx) = mean(Tmax(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 5));
    T.mTMAX_JJA(idx) = mean(Tmax(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 8));
    T.VPD_SON(idx) = VPD(ind-1, 11)- mean(VPD(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 11));
    T.VPD_DJF(idx) = VPD(ind, 2)- mean(VPD(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 2));
    T.VPD_MAM(idx) = VPD(ind, 5)- mean(VPD(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 5));
    T.VPD_JJA(idx) = VPD(ind, 8)- mean(VPD(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 8));
    T.mVPD_SON(idx) = mean(VPD(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 11));
    T.mVPD_DJF(idx) = mean(VPD(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 2));
    T.mVPD_MAM(idx) = mean(VPD(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 5));
    T.mVPD_JJA(idx) = mean(VPD(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 8));
    
    
end
writetable(T, 'StressIndexCovariateTable.csv');

clearvars -except T nmin syear nTopo nSoil nWater nTemp SoilLabs TopoLabs WaterLabs TempLabs AllLabs p smonth;

%% Calibrate LME
% Variables selected from RF model variable importance
FE_TMIN = 'TMIN_DJF';
FE_WB = 'WB_DJF'; 
FE_VPD = 'VPD_JJA';
FE_Site = '(ELEV + TWI + SILT + OC + TN + PHCA + mTMIN_DJF + mWB_DJF + mVPD_JJA)';
REff = ['+ (',FE_TMIN,' | Species) + (',FE_WB,' | Species) + (',FE_VPD,' | Species) + (1 | Species:Site)'];
np = 4; % number of predictors (intercept + main effects)

% Initial full model - selection by BIC
mdl = ['Stress ~ ','(',FE_TMIN,'+',FE_WB,'+',FE_VPD,')*', FE_Site, ' - ', FE_Site];
lme0_bic = fitlme(T, [mdl,REff]);

bic0 = 0;
bic1 = lme0_bic.ModelCriterion.BIC;
lme1 = lme0_bic;

while bic1 < bic0
    
    bic0 = bic1;
    lme0_bic = lme1;
    
    p = lme0_bic.Coefficients.pValue((np+1):end);
    var = lme0_bic.Coefficients.Name((np+1):end);
    
    var2remove = var{p == max(p)};
    mdl = [mdl, ' - ',var2remove];
    
    lme1 = fitlme(T, [mdl,REff]);
    
    bic1 = lme1.ModelCriterion.BIC;
    
end

%% Get random effects
[B,Bnames, stats] = randomEffects(lme0_bic);
stats = dataset2table(stats);
lme0_coeffs = dataset2table(lme0_bic.Coefficients);

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 5];

% Tmin
idx = strcmp(Bnames.Group, 'Species') & strcmp(Bnames.Name, 'TMIN_DJF');
b = lme0_coeffs.Estimate(2);
bse = lme0_coeffs.SE(2);
u = stats.Estimate(idx);
uli = stats.Lower(idx);
uui = stats.Upper(idx);

subplot(3,1,1)
fill([0 34 34 0], [b-bse*1.96 b-bse*1.96 b+bse*1.96 b+bse*1.96], 'k',...
    'EdgeColor','none', 'FaceAlpha',0.1)
hold on;
plot([0 34], [b b], '-', 'LineWidth',1.5, 'Color',[0.4 0.4 0.4])
for i = 1:length(u)
    plot([i i], [b+uli(i) b+uui(i)], 'k-', 'LineWidth',1);
    scatter(i, b+u(i), 10, 'k', 'filled','MarkerFaceColor','w',...
        'LineWidth',1.5, 'MarkerEdgeColor','k')
end
set(gca, 'XLim',[0 34], 'XTick', 1:33, 'TickDir','out', 'XTickLabels','')
box off;
ylabel('\deltaS^{*} / \deltaTMIN_{DJF}');


% CWB
idx = strcmp(Bnames.Group, 'Species') & strcmp(Bnames.Name, 'WB_DJF');
b = lme0_coeffs.Estimate(3);
bse = lme0_coeffs.SE(3);
u = stats.Estimate(idx);
uli = stats.Lower(idx);
uui = stats.Upper(idx);
spc = Bnames.Level(idx);

subplot(3,1,2)
fill([0 34 34 0], [b-bse*1.96 b-bse*1.96 b+bse*1.96 b+bse*1.96], 'k',...
    'EdgeColor','none', 'FaceAlpha',0.1)
hold on;
plot([0 34], [b b], '-', 'LineWidth',1.5, 'Color',[0.4 0.4 0.4])
for i = 1:length(u)
    plot([i i], [b+uli(i) b+uui(i)], 'k-', 'LineWidth',1);
    scatter(i, b+u(i), 10, 'k', 'filled','MarkerFaceColor','w',...
        'LineWidth',1.5, 'MarkerEdgeColor','k')
end
set(gca, 'XLim',[0 34], 'XTick', 1:33, 'TickDir','out', 'XTickLabels','')
box off;
ylabel('\deltaS^{*} / \deltaCWB_{DJF}');


% VPD
idx = strcmp(Bnames.Group, 'Species') & strcmp(Bnames.Name, 'VPD_JJA');
spc = Bnames.Level(idx);
b = lme0_coeffs.Estimate(4);
bse = lme0_coeffs.SE(4);
u = stats.Estimate(idx);
uli = stats.Lower(idx);
uui = stats.Upper(idx);

subplot(3,1,3)
fill([0 34 34 0], [b-bse*1.96 b-bse*1.96 b+bse*1.96 b+bse*1.96], 'k',...
    'EdgeColor','none', 'FaceAlpha',0.1)
hold on;
plot([0 34], [b b], '-', 'LineWidth',1.5, 'Color',[0.4 0.4 0.4])
for i = 1:length(u)
    plot([i i], [b+uli(i) b+uui(i)], 'k-', 'LineWidth',1);
    scatter(i, b+u(i), 10, 'k', 'filled','MarkerFaceColor','w',...
        'LineWidth',1.5, 'MarkerEdgeColor','k')
end
set(gca, 'XLim',[0 34], 'XTick', 1:33, 'TickDir','out', 'XTickLabels',spc)
xtickangle(45)
box off;
ylabel('\deltaS^{*} / \deltaVPD_{JJA}');

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/supplemental-lme-random-effects-species.tif')
close all;

%% Map site-level random effects
load ./data/TREESI_wEnvFactors;
stats = stats(199:end, :);
site1 = stats.Level;
site1 = split(site1, ' ');
site1 = site1(:, 2);
site2 = {TREESI.SITE}';
index = cell2mat(cellfun(@(a) strmatch(a,site2,'exact'),site1,'uniform',false));
lat = [TREESI.LAT]; lat = lat(index);
lon = [TREESI.LON]; lon = lon(index);

states = shaperead('usastatehi','UseGeoCoords', true);
latlim = [24 49];
lonlim = [-125 -67];

clr = wesanderson('fantasticfox1');
grn = make_cmap([1 1 1; clr(3,:); clr(3,:).^3], 5);
prpl = make_cmap([1 1 1; sqrt(clr(4,:)); clr(4,:)], 5);
cbrew = flipud([flipud(grn(2:end, :)); prpl(2:5, :)]);

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 7 4];

ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',[28,34,40,46],'MLineLocation',10,'MeridianLabel','on',...
        'ParallelLabel','on','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Times New Roman', 'GColor',[0.8 0.8 0.8],...
        'FLineWidth',1, 'FontColor',[0.4 0.4 0.4], 'MLabelLocation',20,...
        'MLabelParallel',24.01);
scatterm(lat, lon, 25, stats.Estimate, 'filled', 'MarkerEdgeColor','k');
geoshow(states,'FaceColor','none','LineWidth',0.4)
caxis([-0.4 0.4])
colormap(cbrew)
axis off;
axis image;

cb = colorbar('eastoutside');
cb.TickLength = 0.08;
ylabel(cb, 'Site-level random effect', 'FontSize',11);

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/supplemental-lme-random-effects-site.tif')
close all;


%% Plot model predictions
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 5 4];
yhat = predict(lme0_bic);
clr = wesanderson('fantasticfox1');
clr2 = make_cmap([1 1 1; sqrt(clr(4,:)); clr(4,:)], 11);

hexscatter(T.Stress,yhat, 'res',25, 'drawEdges','true', 'xlim', [0 1], 'ylim',[0 1]);
caxis([0 100])
colormap(clr2(2:end,:));
set(gca, 'XLim', [-0.05 1.05], 'YLim',[-0.05 1.05], 'TickDir','out',...
    'TickLength',[0.02 0.03], 'XTick',0:0.2:1, 'YTick',0:0.2:1);
xlabel('Observed \itS^{*}');
ylabel('Predicted \itS^{*}');
text(0,1, ['R^{2}_{adj} = ', num2str(round(lme0_bic.Rsquared.Adjusted ,2))]);

cb=colorbar('eastoutside');
cb.TickLabels{11} = '100+';
cb.TickLength = 0.07;
ylabel(cb, 'Count', 'FontSize',11);

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/supplemental-lme-model-skill.tif')
close all;


%% Load RELATIVE tree-ring stress estimates and put in table
load ./data/TREESI_wEnvFactors_RSI;

n = sum([TREESI.END] - syear + 1);

T = table('size',[n 2+nTopo+nSoil+nTemp*2+nWater*2+2], 'VariableTypes',...
    {'double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double',...
    'double','double','double','string','string'},...
    'VariableNames',AllLabs);
idx0 = 1;

for i = 1:length(TREESI)
    
    yr = TREESI(i).YEAR;
    ind = yr >= syear;
    idx = idx0:(idx0+sum(ind)-1);
    idx0 = idx0+sum(ind);
    
    T.Stress(idx) = TREESI(i).TREESI(ind);
    T.Year(idx) = yr(ind);
    T.ELEV(idx) = TREESI(i).ELEV;
    T.SLOPE(idx) = TREESI(i).SLOPE;
    T.UAA(idx) = TREESI(i).UAA;
    T.TWI(idx) = TREESI(i).TWI;
    T.CLAY(idx) = TREESI(i).GSDE.CLAY;
    T.OC(idx) = TREESI(i).GSDE.OC;
    T.PHCA(idx) = TREESI(i).GSDE.PHCA;
    T.SAND(idx) = TREESI(i).GSDE.SAND;
    T.SILT(idx) = TREESI(i).GSDE.SILT;
    T.TN(idx) = TREESI(i).GSDE.TN;
    T.VMC2(idx) = TREESI(i).GSDE.VMC2; % volumetric water content at -33 kPa

    T.Site(idx) = TREESI(i).SITE;
    T.Species(idx) = TREESI(i).SPECIES;
    
    [~, ind] = intersect(TREESI(i).PRISM.YEAR, yr(ind));
    
    % 3-month P-PET sums
    P = reshape((TREESI(i).PRISM.PPT - TREESI(i).PRISM.PET)', 1, []);
    windowSize = 3;
    b = ones(1,windowSize);
    a = 1;
    P = filter(b, a, P);
    P = reshape(P, 12, [])';
    T.WB_SON(idx) = P(ind-1, 11)- mean(P(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 11));
    T.WB_DJF(idx) = P(ind, 2)- mean(P(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 2));
    T.WB_MAM(idx) = P(ind, 5)- mean(P(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 5));
    T.WB_JJA(idx) = P(ind, 8)- mean(P(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 8));
    T.mWB_SON(idx) = mean(P(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 11));
    T.mWB_DJF(idx) = mean(P(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 2));
    T.mWB_MAM(idx) = mean(P(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 5));
    T.mWB_JJA(idx) = mean(P(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 8));
    
    % 3-month Tmin, Tmax, and VPD means
    Tmin = reshape(TREESI(i).PRISM.TMIN', 1, []);
    Tmax = reshape(TREESI(i).PRISM.TMAX', 1, []);
    VPD = reshape(TREESI(i).PRISM.VPD', 1, []);
    windowSize = 3;
    b = ones(1,windowSize)/windowSize;
    a = 1;
    Tmin = reshape(filter(b, a, Tmin), 12, [])';
    Tmax = reshape(filter(b, a, Tmax), 12, [])';
    VPD = reshape(filter(b, a, VPD), 12, [])';
    T.TMIN_SON(idx) = Tmin(ind-1, 11)- mean(Tmin(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 11));
    T.TMIN_DJF(idx) = Tmin(ind, 2)- mean(Tmin(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 2));
    T.TMIN_MAM(idx) = Tmin(ind, 5)- mean(Tmin(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 5));
    T.TMIN_JJA(idx) = Tmin(ind, 8)- mean(Tmin(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 8));
    T.mTMIN_SON(idx) = mean(Tmin(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 11));
    T.mTMIN_DJF(idx) = mean(Tmin(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 2));
    T.mTMIN_MAM(idx) = mean(Tmin(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 5));
    T.mTMIN_JJA(idx) = mean(Tmin(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 8));
    T.TMAX_SON(idx) = Tmax(ind-1, 11)- mean(Tmax(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 11));
    T.TMAX_DJF(idx) = Tmax(ind, 2)- mean(Tmax(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 2));
    T.TMAX_MAM(idx) = Tmax(ind, 5)- mean(Tmax(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 5));
    T.TMAX_JJA(idx) = Tmax(ind, 8)- mean(Tmax(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 8));
    T.mTMAX_SON(idx) = mean(Tmax(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 11));
    T.mTMAX_DJF(idx) = mean(Tmax(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 2));
    T.mTMAX_MAM(idx) = mean(Tmax(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 5));
    T.mTMAX_JJA(idx) = mean(Tmax(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 8));
    T.VPD_SON(idx) = VPD(ind-1, 11)- mean(VPD(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 11));
    T.VPD_DJF(idx) = VPD(ind, 2)- mean(VPD(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 2));
    T.VPD_MAM(idx) = VPD(ind, 5)- mean(VPD(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 5));
    T.VPD_JJA(idx) = VPD(ind, 8)- mean(VPD(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 8));
    T.mVPD_SON(idx) = mean(VPD(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 11));
    T.mVPD_DJF(idx) = mean(VPD(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 2));
    T.mVPD_MAM(idx) = mean(VPD(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 5));
    T.mVPD_JJA(idx) = mean(VPD(TREESI(i).PRISM.YEAR>=1981 & TREESI(i).PRISM.YEAR<=2010, 8));
    
    
end
writetable(T, 'StressIndexCovariateTable_RSI.csv');

clearvars -except T nmin syear nTopo nSoil nWater nTemp SoilLabs TopoLabs WaterLabs TempLabs AllLabs p smonth;

%% Calibrate LME
% Variables selected from RF model variable importance
FE_TMIN = 'TMIN_DJF';
FE_WB = 'WB_DJF'; 
FE_VPD = 'VPD_JJA';
FE_Site = '(ELEV + TWI + SILT + OC + TN + PHCA + mTMIN_DJF + mWB_DJF + mVPD_JJA)';
REff = ['+ (',FE_TMIN,' | Species) + (',FE_WB,' | Species) + (',FE_VPD,' | Species) + (1 | Species:Site)'];
np = 4; % number of predictors (intercept + main effects)

% Initial full model - selection by BIC
mdl = ['Stress ~ ','(',FE_TMIN,'+',FE_WB,'+',FE_VPD,')*', FE_Site, ' - ', FE_Site];
lme0_bic = fitlme(T, [mdl,REff]);

bic0 = 0;
bic1 = lme0_bic.ModelCriterion.BIC;
lme1 = lme0_bic;

while bic1 < bic0
    
    bic0 = bic1;
    lme0_bic = lme1;
    
    p = lme0_bic.Coefficients.pValue((np+1):end);
    var = lme0_bic.Coefficients.Name((np+1):end);
    
    var2remove = var{p == max(p)};
    mdl = [mdl, ' - ',var2remove];
    
    lme1 = fitlme(T, [mdl,REff]);
    
    bic1 = lme1.ModelCriterion.BIC;
    
end

coeffs = dataset2table(lme0_bic.Coefficients);

%% Plot model predictions
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 5 4];
yhat = predict(lme0_bic);
clr = wesanderson('fantasticfox1');
clr2 = make_cmap([1 1 1; sqrt(clr(4,:)); clr(4,:)], 11);

hexscatter(T.Stress,yhat, 'res',50, 'drawEdges','true', 'xlim', [0 2], 'ylim',[0 2]);
caxis([0 100])
colormap(clr2(2:end,:));
set(gca, 'XLim', [-0.1 2.1], 'YLim',[-0.1 2.1], 'TickDir','out',...
    'TickLength',[0.02 0.03], 'XTick',0:0.5:2, 'YTick',0:0.5:2);
xlabel('Observed \itS_{r}');
ylabel('Predicted \itS_{r}');
text(0,2, ['R^{2}_{adj} = ', num2str(round(lme0_bic.Rsquared.Adjusted ,2))]);

cb=colorbar('eastoutside');
cb.TickLabels{11} = '100+';
cb.TickLength = 0.07;
ylabel(cb, 'Count', 'FontSize',11);

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/supplemental-lme-model-skill-rsi.tif')
close all;


