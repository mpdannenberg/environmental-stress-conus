% Make variable importance plot from RF models of relative stress index

n_min = 75;
syear = 1971;
clr = wesanderson('fantasticfox1');
clr2 = make_cmap([1 1 1; clr(4,:)], 8);

% Number of vars for each environmental stress type
nTopo = 4;
nSoil = 6; % Subject to change
nTemp = 8;
nWater = 8;

% Names of variables
TopoLabs = {'ELEV','SLOPE','UAA','TWI'};
SoilLabs = {'CLAY','OC','PHCA','SAND','SILT','TN'};
TempLabs = {'TMIN_{SON}','TMIN_{DJF}','TMIN_{MAM}','TMIN_{JJA}',...
	'TMAX_{SON}','TMAX_{DJF}','TMAX_{MAM}','TMAX_{JJA}'};
WaterLabs = {'WB_{SON}','WB_{DJF}','WB_{MAM}','WB_{JJA}',...
	'VPD_{SON}','VPD_{DJF}','VPD_{MAM}','VPD_{JJA}'};
AllLabs = [TopoLabs SoilLabs TempLabs WaterLabs];

p = length(AllLabs);

%% Load/process datasets
load ./data/TREESI_wEnvFactors_RSI;
load ./data/RF_v6_RSI_TSC.mat;

EcoL1 = [TREESI.EcoL1];
EcoL2 = [TREESI.EcoL2];
EcoL3 = [TREESI.EcoL3];

nyrs = [TREESI.END] - (syear-1);

L1list = sort(unique(EcoL1));
L1list = L1list(L1list >0);
for i = 1:length(L1list)
    L1n(i) = sum(nyrs(EcoL1 == L1list(i)));
end
L1list = L1list(L1n >= n_min); L1n = L1n(L1n >= n_min);

L2list = sort(unique(EcoL2));
L2list = L2list(L2list >0);
for i = 1:length(L2list)
    L2n(i) = sum(nyrs(EcoL2 == L2list(i)));
end
L2list = L2list(L2n >= n_min); L2n = L2n(L2n >= n_min);

L3list = sort(unique(EcoL3));
L3list = L3list(L3list >0);
for i = 1:length(L3list)
    L3n(i) = sum(nyrs(EcoL3 == L3list(i)));
end
L3list = L3list(L3n >= n_min); L3n = L3n(L3n >= n_min);

EcoL1Labs = {'Northern Forests','Northwestern Forested Mountains',...
    'Marine West Coast Forest','Eastern Temperate Forests','Great Plains',...
    'North American Deserts','Mediterranean California',...
    'Southern Semi-Arid Highlands','Temperate Sierras'};
EcoL2Labs = cellstr(num2str(L2list'/10));
EcoL3Labs = cellstr(strcat(num2str( floor(L3list'/100)/10 ),repmat('.',length(L3list),1), num2str( 100*(L3list'/100-floor(L3list'/100)))));
EcoL3Labs = regexprep(EcoL3Labs, ' ', ''); % remove extra spaces

%% Get variable importance for each L1 model
VarImp = NaN(length(AllLabs), length(L1list));
for i = 1:length(L1list)
    
    eval(['mdl = B_', num2str(L1list(i)),';']);
    VarImp(:, i) = mdl.OOBPermutedPredictorDeltaError;
    
end
wgtMean = sum(VarImp .* repmat(L1n/sum(L1n), length(AllLabs), 1), 2);
VarImp(:, end+1) = wgtMean;

% Show variable importance
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 2];

imagesc(1:length(AllLabs), 1:(length(L1list)+1), VarImp');
caxis([0 2])
colormap(clr2)
set(gca, 'XTick', 1:length(AllLabs), 'XTickLabels',AllLabs,...
    'YTick',1:(length(L1list)+1), 'YTickLabels',[EcoL1Labs {'\bfWeighted Mean'}], 'FontSize',7,...
    'TickDir','out')
hold on;
for i = 1:length(AllLabs)
    plot([i+0.5 i+0.5], [0 length(L1list)+2], 'k-')
end
for i = 0:length(L1list)
    plot([0 length(AllLabs)+2], [i+0.5 i+0.5], 'k-')
end
box off;
xtickangle(360-45)
ax = gca;
ax.Position = [0.2369 0.2565 0.64 0.6685];

cb = colorbar('eastoutside');
cb.Position = [0.89 0.2552 0.0342 0.6719];
ylb = ylabel(cb, '\DeltaError');
cb.FontSize = 9;
cb.Ticks = 0:0.25:2;
cb.TickLabels = {'0','','0.5','','1','','1.5','','2'};
cb.TickLength = 0.17;
ylb.Position(1) = 2.2;

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/supplemental-rsi-rf-variable-importance-L1.tif')
close all;

%% Get variable importance for each L2 model
VarImp = NaN(length(AllLabs), length(L2list));
for i = 1:length(L2list)
    
    eval(['mdl = B_', num2str(L2list(i)),';']);
    VarImp(:, i) = mdl.OOBPermutedPredictorDeltaError;
    
end
wgtMean = sum(VarImp .* repmat(L2n/sum(L2n), length(AllLabs), 1), 2);
VarImp(:, end+1) = wgtMean;

% Show variable importance
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 3.25];

imagesc(1:length(AllLabs), 1:(length(L2list)+1), VarImp');
caxis([0 2])
colormap(clr2)
set(gca, 'XTick', 1:length(AllLabs), 'XTickLabels',AllLabs,...
    'YTick',1:(length(L2list)+1), 'YTickLabels',[EcoL2Labs' {'\bfWeighted Mean'}], 'FontSize',7,...
    'TickDir','out')
hold on;
for i = 1:length(AllLabs)
    plot([i+0.5 i+0.5], [0 length(L2list)+2], 'k-')
end
for i = 0:length(L2list)
    plot([0 length(AllLabs)+2], [i+0.5 i+0.5], 'k-')
end
box off;
xtickangle(360-45)
ax = gca;
ax.Position = [0.1300    0.1585    0.74    0.7665];

cb = colorbar('eastoutside');
cb.Position = [0.89    0.1571    0.0342    0.7692];
ylb = ylabel(cb, '\DeltaError');
cb.FontSize = 9;
cb.Ticks = 0:0.25:2;
cb.TickLabels = {'0','','0.5','','1','','1.5','','2'};
cb.TickLength = 0.09;
ylb.Position(1) = 2.1;

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/supplemental-rsi-rf-variable-importance-L2.tif')
close all;

%% Get variable importance for each L3 model
VarImp = NaN(length(AllLabs), length(L3list));
for i = 1:length(L3list)
    
    eval(['mdl = B_', num2str(L3list(i)),';']);
    VarImp(:, i) = mdl.OOBPermutedPredictorDeltaError;
    
end
wgtMean = sum(VarImp .* repmat(L3n/sum(L3n), length(AllLabs), 1), 2);
VarImp(:, end+1) = wgtMean;

% Show variable importance
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 0 6.5 8.5];

imagesc(1:length(AllLabs), 1:(length(L3list)+1), VarImp');
caxis([0 2])
colormap(clr2)
set(gca, 'XTick', 1:length(AllLabs), 'XTickLabels',AllLabs,...
    'YTick',1:(length(L3list)+1), 'YTickLabels',[EcoL3Labs' {'\bfWeighted Mean'}], 'FontSize',7,...
    'TickDir','out')
hold on;
for i = 1:length(AllLabs)
    plot([i+0.5 i+0.5], [0 length(L3list)+2], 'k-')
end
for i = 0:length(L3list)
    plot([0 length(AllLabs)+2], [i+0.5 i+0.5], 'k-')
end
box off;
xtickangle(360-45)
ax = gca;
ax.Position = [0.1300    0.07    0.7750    0.84];

cb = colorbar('northoutside');
cb.Position = [0.1298    0.925    0.7756    0.0261];
ylb = ylabel(cb, '\DeltaError');
cb.FontSize = 9;
cb.Ticks = 0:0.25:2;
cb.TickLabels = {'0','','0.5','','1','','1.5','','2'};
cb.TickLength = 0.045;

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/supplemental-rsi-rf-variable-importance-L3.tif')
close all;
