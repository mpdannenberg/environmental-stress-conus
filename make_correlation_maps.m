% Map spearman rank correlations between climate variables and
% environmental stress index

alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

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

load ./data/TREESI_wEnvFactors;

EcoL1 = [TREESI.EcoL1];
EcoL2 = [TREESI.EcoL2];
EcoL3 = [TREESI.EcoL3];

L1list = sort(unique(EcoL1));
L1list = L1list(L1list >0);

L2list = sort(unique(EcoL2));
L2list = L2list(L2list >0);

L3list = sort(unique(EcoL3));
L3list = L3list(L3list >0);

load ./data/RF_v4_TSC_stats;
[EcoL1, R] = geotiffread('./data/us_EcoL1_4km.tif');
[EcoL2] = geotiffread('./data/us_EcoL2_4km.tif');
[EcoL3] = geotiffread('./data/us_EcoL3_4km.tif');
load ./data/prism_latlon;
load ./data/conus_mask;

states = shaperead('usastatehi','UseGeoCoords', true);
latlim = [25 50];
lonlim = [-126 -65];

clr = wesanderson('fantasticfox1');
grn = make_cmap([1 1 1; clr(3,:); clr(3,:).^3], 5);
prpl = make_cmap([1 1 1; sqrt(clr(4,:)); clr(4,:)], 5);
cbrew = flipud([flipud(grn(2:end, :)); prpl(2:5, :)]);

% Temperature
h=figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 5 6.33];
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
    
    r(p>0.05) = NaN;
    
    subplot(4,2,pos(i))
    if i~=5
        ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
            'on','PLineLocation',10,'MLineLocation',20,'MeridianLabel','off',...
            'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
            'none', 'FontName', 'Helvetica',...
            'FLineWidth',1, 'FontSize',7, 'MLabelParallel',25.01,...
            'GColor',[0.8 0.8 0.8], 'FontColor',[0.4 0.4 0.4]);
    else
        ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
            'on','PLineLocation',10,'MLineLocation',20,'MeridianLabel','on',...
            'ParallelLabel','on','GLineWidth',0.3,'Frame','off','FFaceColor',...
            'none', 'FontName', 'Helvetica',...
            'FLineWidth',1, 'FontSize',7, 'MLabelParallel',25.01,...
            'GColor',[0.8 0.8 0.8], 'FontColor',[0.4 0.4 0.4], 'PLabelMeridian','east');
    end
    surfm(PRISMlat, PRISMlon, r.*CONUS);
    caxis([-0.8 0.8]);
    colormap(cbrew);
    geoshow(states,'FaceColor','none','LineWidth',0.3, 'EdgeColor',[0.6 0.6 0.6])
    axis off;
    axis image;
    if i ==1
        ttl = title('TMIN', 'FontName','Helvetica', 'FontSize',12);
        ttl.Position(2) = ttl.Position(2)+0.05;
    elseif i ==5
        ttl = title('TMAX', 'FontName','Helvetica', 'FontSize',12);
        ttl.Position(2) = ttl.Position(2)+0.05;
    end
    subplotsqueeze(gca, 1.2)
    if i==1
        text(-0.5,0.65,'SON','FontSize',12, 'Rotation',90, 'FontWeight','bold');
    elseif i==2
        text(-0.5,0.65,'DJF','FontSize',12, 'Rotation',90, 'FontWeight','bold');
    elseif i==3
        text(-0.5,0.65,'MAM','FontSize',12, 'Rotation',90, 'FontWeight','bold');
    elseif i==4
        text(-0.5,0.65,'JJA','FontSize',12, 'Rotation',90, 'FontWeight','bold');
    end
    text(-0.38,0.85,alphabet(pos(i)),'FontSize',12, 'FontWeight','bold');
        
end
h = colorbar('southoutside');
h.Position = [0.15 0.06 0.7 0.02];
h.FontName = 'Helvetica';
h.Ticks = -0.8:0.2:0.8;
h.TickLength = 0.038;
h.Label.String = 'Spearman''s \rho';

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/correlations-temperature.tif')
close all;


% Water
h=figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 5 6.33];
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
    
    r(p>0.05) = NaN;
    
    subplot(4,2,pos(i))
    if i~=5
        ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
            'on','PLineLocation',10,'MLineLocation',20,'MeridianLabel','off',...
            'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
            'none', 'FontName', 'Helvetica',...
            'FLineWidth',1, 'FontSize',7, 'MLabelParallel',25.01,...
            'GColor',[0.8 0.8 0.8], 'FontColor',[0.4 0.4 0.4]);
    else
        ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
            'on','PLineLocation',10,'MLineLocation',20,'MeridianLabel','on',...
            'ParallelLabel','on','GLineWidth',0.3,'Frame','off','FFaceColor',...
            'none', 'FontName', 'Helvetica',...
            'FLineWidth',1, 'FontSize',7, 'MLabelParallel',25.01,...
            'GColor',[0.8 0.8 0.8], 'FontColor',[0.4 0.4 0.4], 'PLabelMeridian','east');
    end
    surfm(PRISMlat, PRISMlon, r.*CONUS);
    caxis([-0.8 0.8]);
    colormap(cbrew);
    geoshow(states,'FaceColor','none','LineWidth',0.3, 'EdgeColor',[0.6 0.6 0.6])
    axis off;
    axis image;
    if i ==1
        ttl = title('CWB', 'FontName','Helvetica', 'FontSize',12);
        ttl.Position(2) = ttl.Position(2)+0.05;
    elseif i ==5
        ttl = title('VPD', 'FontName','Helvetica', 'FontSize',12);
        ttl.Position(2) = ttl.Position(2)+0.05;
    end
    subplotsqueeze(gca, 1.2)
    if i==1
        text(-0.5,0.65,'SON','FontSize',12, 'Rotation',90, 'FontWeight','bold');
    elseif i==2
        text(-0.5,0.65,'DJF','FontSize',12, 'Rotation',90, 'FontWeight','bold');
    elseif i==3
        text(-0.5,0.65,'MAM','FontSize',12, 'Rotation',90, 'FontWeight','bold');
    elseif i==4
        text(-0.5,0.65,'JJA','FontSize',12, 'Rotation',90, 'FontWeight','bold');
    end
    text(-0.38,0.85,alphabet(pos(i)),'FontSize',12, 'FontWeight','bold');
        
end
h = colorbar('southoutside');
h.Position = [0.15 0.06 0.7 0.02];
h.FontName = 'Helvetica';
h.Ticks = -0.8:0.2:0.8;
h.TickLength = 0.038;
h.Label.String = 'Spearman''s \rho';

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/correlations-water.tif')
close all;
