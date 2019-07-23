% Make variable importance plot from RF models

load ./data/RF_v4_TSC.mat;

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
load ./data/TREESI_wEnvFactors;

EcoL1 = [TREESI.EcoL1];
EcoL2 = [TREESI.EcoL2];
EcoL3 = [TREESI.EcoL3];

L1list = sort(unique(EcoL1));
L1list = L1list(L1list >0);
for i = 1:length(L1list)
    L1n(i) = sum(EcoL1 == L1list(i));
end

L2list = sort(unique(EcoL2));
L2list = L2list(L2list >0);

L3list = sort(unique(EcoL3));
L3list = L3list(L3list >0);

EcoLabs = {'Northern Forests','Northwestern Forested Mountains',...
    'Marine West Coast Forest','Eastern Temperate Forests','Great Plains',...
    'North American Deserts','Mediterranean California',...
    'Southern Semi-Arid Highlands','Temperate Sierras'};

%% Get variable importance from each model
VarImp = NaN(length(AllLabs), length(L1list));
for i = 1:length(L1list)
    
    eval(['mdl = B_', num2str(L1list(i)),';']);
    VarImp(:, i) = mdl.OOBPermutedPredictorDeltaError;
    
end
wgtMean = sum(VarImp .* repmat(L1n/sum(L1n), length(AllLabs), 1), 2);
VarImp(:, end+1) = wgtMean;

%% Show variable importance
clr = wesanderson('fantasticfox1');
clr2 = make_cmap([1 1 1; clr(4,:)], 9);

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 2];

imagesc(1:length(AllLabs), 1:(length(L1list)+1), VarImp');
caxis([0 1.25])
colormap(clr2)
set(gca, 'XTick', 1:length(AllLabs), 'XTickLabels',AllLabs,...
    'YTick',1:(length(L1list)+1), 'YTickLabels',[EcoLabs {'\bfWeighted Mean'}], 'FontSize',7,...
    'TickDir','out')
xtickangle(360-45)
ax = gca;
ax.Position = [0.2369 0.2565 0.64 0.6685];

cb = colorbar('eastoutside');
cb.Position = [0.9 0.2552 0.0342 0.6719];
ylb = ylabel(cb, '\DeltaError');
cb.FontSize = 9;
ylb.Position(1) = 2.1;

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/rf-variable-importance.tif')
close all;
