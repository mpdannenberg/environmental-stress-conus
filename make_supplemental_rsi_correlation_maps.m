%%
%% NEED TO UPDATE THIS FOR CORRELATIONS INSTEAD OF IMPORTANCE
%%

% Map variable importance

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

load TREESI_wEnvFactors3;

EcoL1 = [TREESI.EcoL1];
EcoL2 = [TREESI.EcoL2];
EcoL3 = [TREESI.EcoL3];

L1list = sort(unique(EcoL1));
L1list = L1list(L1list >0);

L2list = sort(unique(EcoL2));
L2list = L2list(L2list >0);

L3list = sort(unique(EcoL3));
L3list = L3list(L3list >0);


cd RF_v6_RSI_TSC;
load RF_v6_RSI_TSC_stats;
cd ..;

cd EcoRegions;
[EcoL1, R] = geotiffread('us_EcoL1_4km.tif');
[EcoL2] = geotiffread('us_EcoL2_4km.tif');
[EcoL3] = geotiffread('us_EcoL3_4km.tif');
cd ../PRISM;
load prism_latlon;
cd ..;
load conus_mask;

states = shaperead('usastatehi','UseGeoCoords', true);
latlim = [25 50];
lonlim = [-126 -65];

cbrew = [140,81,10
    191,129,45
    223,194,125
    246,232,195
    199,234,229
    128,205,193
    53,151,143
    1,102,94]/255;


% Temperature
ClimLabs = [TempLabs WaterLabs];
h=figure('Color','w');
h.Units = 'inches';
h.Position = [2 3 5 6.33];
pos = [1:2:7 2:2:8];
for i=1:8
    
    r = NaN(size(EcoL2));
    p = NaN(size(EcoL2));
    
    % Add L1 stats
    for j=1:length(L1list)
        el = L1list(j);
        rel = L1_corrs(j, i+nTopo+nSoil, 1);
        pel = L1_corrs(j, i+nTopo+nSoil, 2);
        
        if ~isnan(rel)
            r(EcoL1==el) = rel;
            p(EcoL1==el) = pel;
        end
    end
    
    % Add L2 stats
    for j=1:length(L2list)
        el = L2list(j);
        rel = L2_corrs(j, i+nTopo+nSoil, 1);
        pel = L2_corrs(j, i+nTopo+nSoil, 2);
        
        if ~isnan(rel)
            r(EcoL2==el) = rel;
            p(EcoL2==el) = pel;
        end
    end
    
    % Add L3 stats
    for j=1:length(L3list)
        el = L3list(j);
        rel = L3_corrs(j, i+nTopo+nSoil, 1);
        pel = L3_corrs(j, i+nTopo+nSoil, 2);
        
        if ~isnan(rel)
            r(EcoL3==el) = rel;
            p(EcoL3==el) = pel;
        end
    end
    
    %p2 = double(p<0.05).*CONUS;
    %p2(p2==0) = NaN;
    r(p>0.05) = NaN;
    
    subplot(4,2,pos(i))
    % map...
    if i~=1
        ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
            'on','PLineLocation',10,'MLineLocation',20,'MeridianLabel','off',...
            'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
            'none', 'FontName', 'Times New Roman',...
            'FLineWidth',1, 'FontSize',7, 'MLabelParallel',25.01,...
            'GColor',[0.8 0.8 0.8], 'FontColor',[0.4 0.4 0.4]);
    else
        ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
            'on','PLineLocation',10,'MLineLocation',20,'MeridianLabel','on',...
            'ParallelLabel','on','GLineWidth',0.3,'Frame','off','FFaceColor',...
            'none', 'FontName', 'Times New Roman',...
            'FLineWidth',1, 'FontSize',7, 'MLabelParallel',25.01,...
            'GColor',[0.8 0.8 0.8], 'FontColor',[0.4 0.4 0.4]);
    end
    surfm(PRISMlat, PRISMlon, r.*CONUS);
    caxis([-0.8 0.8]);
    colormap(cbrew);
    geoshow(states,'FaceColor','none','LineWidth',0.3, 'EdgeColor',[0.6 0.6 0.6])
    %contourm(PRISMlat, PRISMlon, p2, 'LineColor','k');
%     plotm(LAT, LON, 'k^', 'MarkerSize',4, 'MarkerFaceColor','w')
%     plotm(tr_lat, tr_lon, 'k^', 'MarkerSize',4, 'MarkerFaceColor','k')
    axis off;
    axis image;
    title(TempLabs{i}, 'FontName','Times New Roman');
end
h = colorbar('southoutside');
h.Position = [0.417 0.07 0.2 0.02];
h.FontName = 'Times New Roman';
h.Ticks = -0.8:0.4:0.8;
h.Label.String = 'Spearman''s \rho';

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r600','RF_v6_RSI_TSC_TempCorrs.tif')
close all;


% Water
ClimLabs = [TempLabs WaterLabs];
h=figure('Color','w');
h.Units = 'inches';
h.Position = [2 3 5 6.33];
pos = [1:2:7 2:2:8];
for i=1:8
    
    r = NaN(size(EcoL2));
    p = NaN(size(EcoL2));
    
    % Add L1 stats
    for j=1:length(L1list)
        el = L1list(j);
        rel = L1_corrs(j, i+nTopo+nSoil+nTemp, 1);
        pel = L1_corrs(j, i+nTopo+nSoil+nTemp, 2);
        
        if ~isnan(rel)
            r(EcoL1==el) = rel;
            p(EcoL1==el) = pel;
        end
    end
    
    % Add L2 stats
    for j=1:length(L2list)
        el = L2list(j);
        rel = L2_corrs(j, i+nTopo+nSoil+nTemp, 1);
        pel = L2_corrs(j, i+nTopo+nSoil+nTemp, 2);
        
        if ~isnan(rel)
            r(EcoL2==el) = rel;
            p(EcoL2==el) = pel;
        end
    end
    
    % Add L3 stats
    for j=1:length(L3list)
        el = L3list(j);
        rel = L3_corrs(j, i+nTopo+nSoil+nTemp, 1);
        pel = L3_corrs(j, i+nTopo+nSoil+nTemp, 2);
        
        if ~isnan(rel)
            r(EcoL3==el) = rel;
            p(EcoL3==el) = pel;
        end
    end
    
    %p2 = double(p<0.05).*CONUS;
    %p2(p2==0) = NaN;
    r(p>0.05) = NaN;
    
    subplot(4,2,pos(i))
    % map...
    if i~=1
        ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
            'on','PLineLocation',10,'MLineLocation',20,'MeridianLabel','off',...
            'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
            'none', 'FontName', 'Times New Roman',...
            'FLineWidth',1, 'FontSize',7, 'MLabelParallel',25.01,...
            'GColor',[0.8 0.8 0.8], 'FontColor',[0.4 0.4 0.4]);
    else
        ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
            'on','PLineLocation',10,'MLineLocation',20,'MeridianLabel','on',...
            'ParallelLabel','on','GLineWidth',0.3,'Frame','off','FFaceColor',...
            'none', 'FontName', 'Times New Roman',...
            'FLineWidth',1, 'FontSize',7, 'MLabelParallel',25.01,...
            'GColor',[0.8 0.8 0.8], 'FontColor',[0.4 0.4 0.4]);
    end
    surfm(PRISMlat, PRISMlon, r.*CONUS);
    caxis([-0.8 0.8]);
    colormap(cbrew);
    geoshow(states,'FaceColor','none','LineWidth',0.3, 'EdgeColor',[0.6 0.6 0.6])
    %contourm(PRISMlat, PRISMlon, p2, 'LineColor','k');
%     plotm(LAT, LON, 'k^', 'MarkerSize',4, 'MarkerFaceColor','w')
%     plotm(tr_lat, tr_lon, 'k^', 'MarkerSize',4, 'MarkerFaceColor','k')
    axis off;
    axis image;
    title(WaterLabs{i}, 'FontName','Times New Roman');
end
h = colorbar('southoutside');
h.Position = [0.417 0.07 0.2 0.02];
h.FontName = 'Times New Roman';
h.Ticks = -0.8:0.4:0.8;
h.Label.String = 'Spearman''s \rho';

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r600','RF_v6_RSI_TSC_WaterCorrs.tif')
close all;
% 
% 
% 
% % Climate
% ClimLabs = [TempLabs WaterLabs];
% h=figure('Color','w');
% h.Units = 'inches';
% h.Position = [2 4 7.3 6];
% for i=1:16
%     
%     r = NaN(size(EcoL2));
%     p = NaN(size(EcoL2));
%     
%     % Add L1 stats
%     for j=1:length(L1list)
%         el = L1list(j);
%         rel = L1_corrs(j, i+nTopo+nSoil, 1);
%         pel = L1_corrs(j, i+nTopo+nSoil, 2);
%         
%         if ~isnan(rel)
%             r(EcoL1==el) = rel;
%             p(EcoL1==el) = pel;
%         end
%     end
%     
%     % Add L2 stats
%     for j=1:length(L2list)
%         el = L2list(j);
%         rel = L2_corrs(j, i+nTopo+nSoil, 1);
%         pel = L2_corrs(j, i+nTopo+nSoil, 2);
%         
%         if ~isnan(rel)
%             r(EcoL2==el) = rel;
%             p(EcoL2==el) = pel;
%         end
%     end
%     
%     % Add L3 stats
%     for j=1:length(L3list)
%         el = L3list(j);
%         rel = L3_corrs(j, i+nTopo+nSoil, 1);
%         pel = L3_corrs(j, i+nTopo+nSoil, 2);
%         
%         if ~isnan(rel)
%             r(EcoL3==el) = rel;
%             p(EcoL3==el) = pel;
%         end
%     end
%     
%     %p2 = double(p<0.05).*CONUS;
%     %p2(p2==0) = NaN;
%     r(p>0.05) = NaN;
%     
%     subplot(4,4,i)
%     % map...
%     if i~=13
%         ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
%             'on','PLineLocation',5,'MLineLocation',10,'MeridianLabel','off',...
%             'ParallelLabel','off','GLineWidth',0.5,'Frame','off','FFaceColor',...
%             'none', 'FontName', 'Times New Roman',...
%             'FLineWidth',1);
%     else
%         ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
%             'on','PLineLocation',5,'MLineLocation',10,'MeridianLabel','on',...
%             'ParallelLabel','on','GLineWidth',0.5,'Frame','off','FFaceColor',...
%             'none', 'FontName', 'Times New Roman',...
%             'FLineWidth',1);
%     end
%     surfm(PRISMlat, PRISMlon, r.*CONUS);
%     caxis([-0.8 0.8]);
%     colormap(cbrew);
%     geoshow(states,'FaceColor','none','LineWidth',0.3, 'EdgeColor',[0.6 0.6 0.6])
%     %contourm(PRISMlat, PRISMlon, p2, 'LineColor','k');
% %     plotm(LAT, LON, 'k^', 'MarkerSize',4, 'MarkerFaceColor','w')
% %     plotm(tr_lat, tr_lon, 'k^', 'MarkerSize',4, 'MarkerFaceColor','k')
%     axis off;
%     axis image;
%     title(ClimLabs{i}, 'FontName','Times New Roman');
% end
% h = colorbar('southoutside');
% h.Position = [0.417 0.33 0.2 0.05];
% h.FontName = 'Times New Roman';
% h.Ticks = -0.8:0.4:0.8;
% h.Label.String = 'Spearman''s \rho';
% 
% set(gcf,'PaperPositionMode','auto')
% print('-dtiff','-f1','-r600','RF_v6_RSI_TSC_ClimCorrs.tif')
% close all;
% 




