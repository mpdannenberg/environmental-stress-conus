% Make maps comparing TSC model skill to C model skill

[EcoL1, R] = geotiffread('./data/us_EcoL1_4km.tif');
[EcoL2] = geotiffread('./data/us_EcoL2_4km.tif');
[EcoL3] = geotiffread('./data/us_EcoL3_4km.tif');
load ./data/prism_latlon;
load ./data/conus_mask;

states = shaperead('usastatehi','UseGeoCoords', true);
latlim = [25 50];
lonlim = [-126 -65];

clr = wesanderson('fantasticfox1');
cbrew = make_cmap([1 1 1; clr(4,:)], 8);

h=figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 4];

%% Absolute Stress Index (TSC)
load ./data/RF_v4_TSC_stats;

r = NaN(size(EcoL2));

% Add L1 stats
n = size(L1_stats, 1);
for j=1:n
    el = L1_stats(j, 1);
    rel = L1_stats(j, 3);

    if ~isnan(rel)
        r(EcoL1==el) = rel;
    end
end

% Add L2 stats
n = size(L2_stats, 1);
for j=1:n
    el = L2_stats(j, 1);
    rel = L2_stats(j, 3);

    if ~isnan(rel)
        r(EcoL2==el) = rel;
    end
end

% Add L3 stats
n = size(L3_stats, 1);
for j=1:n
    el = L3_stats(j, 1);
    rel = L3_stats(j, 3);

    if ~isnan(rel)
        r(EcoL3==el) = rel;
    end
end

subplot(2,2,1)
ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
    'on','PLineLocation',10,'MLineLocation',10,'MeridianLabel','on',...
    'ParallelLabel','on','GLineWidth',0.3,'Frame','off','FFaceColor',...
    'none', 'FontName', 'Times New Roman',...
    'FLineWidth',1, 'FontSize',7, 'MLabelParallel',25.01,...
    'GColor',[0.8 0.8 0.8], 'FontColor',[0.4 0.4 0.4], 'PLabelMeridian','east');
surfm(PRISMlat, PRISMlon, r.*CONUS);
caxis([0 1]);
colormap(cbrew);
geoshow(states,'FaceColor','none','LineWidth',0.3, 'EdgeColor',[0.4 0.4 0.4])
axis off;
axis image;
text(-0.42, 0.85, 'A', 'FontName','Times New Roman', 'FontSize',12, 'FontWeight','bold');
text(-0.05, 0.96, 'TSC model', 'FontName','Times New Roman', 'FontSize',14,...
    'FontWeight','bold','HorizontalAlignment','center');
ax.Position = [0.05 0.56 0.4 0.36];

%% Absolute Stress Index (C)
load ./data/RF_v4_C_stats;

r = NaN(size(EcoL2));

% Add L1 stats
n = size(L1_stats, 1);
for j=1:n
    el = L1_stats(j, 1);
    rel = L1_stats(j, 3);

    if ~isnan(rel)
        r(EcoL1==el) = rel;
    end
end

% Add L2 stats
n = size(L2_stats, 1);
for j=1:n
    el = L2_stats(j, 1);
    rel = L2_stats(j, 3);

    if ~isnan(rel)
        r(EcoL2==el) = rel;
    end
end

% Add L3 stats
n = size(L3_stats, 1);
for j=1:n
    el = L3_stats(j, 1);
    rel = L3_stats(j, 3);

    if ~isnan(rel)
        r(EcoL3==el) = rel;
    end
end

subplot(2,2,2)
ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
    'on','PLineLocation',10,'MLineLocation',10,'MeridianLabel','off',...
    'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
    'none', 'FontName', 'Times New Roman',...
    'FLineWidth',1, 'FontSize',7, 'MLabelParallel',25.01,...
    'GColor',[0.8 0.8 0.8], 'FontColor',[0.4 0.4 0.4]);
surfm(PRISMlat, PRISMlon, r.*CONUS);
caxis([0 1]);
colormap(cbrew);
geoshow(states,'FaceColor','none','LineWidth',0.3, 'EdgeColor',[0.4 0.4 0.4])
axis off;
axis image;
text(-0.42, 0.85, 'B', 'FontName','Times New Roman', 'FontSize',12, 'FontWeight','bold');
text(-0.05, 0.96, 'C model', 'FontName','Times New Roman', 'FontSize',14,...
    'FontWeight','bold','HorizontalAlignment','center');
ax.Position = [0.5703 0.56 0.4 0.36];


%% Relative Stress Index (TSC)
load ./data/RF_v6_RSI_TSC_stats;

r = NaN(size(EcoL2));

% Add L1 stats
n = size(L1_stats, 1);
for j=1:n
    el = L1_stats(j, 1);
    rel = L1_stats(j, 3);

    if ~isnan(rel)
        r(EcoL1==el) = rel;
    end
end

% Add L2 stats
n = size(L2_stats, 1);
for j=1:n
    el = L2_stats(j, 1);
    rel = L2_stats(j, 3);

    if ~isnan(rel)
        r(EcoL2==el) = rel;
    end
end

% Add L3 stats
n = size(L3_stats, 1);
for j=1:n
    el = L3_stats(j, 1);
    rel = L3_stats(j, 3);

    if ~isnan(rel)
        r(EcoL3==el) = rel;
    end
end

subplot(2,2,3)
ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
    'on','PLineLocation',10,'MLineLocation',10,'MeridianLabel','off',...
    'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
    'none', 'FontName', 'Times New Roman',...
    'FLineWidth',1, 'FontSize',7, 'MLabelParallel',25.01,...
    'GColor',[0.8 0.8 0.8], 'FontColor',[0.4 0.4 0.4]);
surfm(PRISMlat, PRISMlon, r.*CONUS);
caxis([0 1]);
colormap(cbrew);
geoshow(states,'FaceColor','none','LineWidth',0.3, 'EdgeColor',[0.4 0.4 0.4])
axis off;
axis image;
text(-0.42, 0.85, 'C', 'FontName','Times New Roman', 'FontSize',12, 'FontWeight','bold');
ax.Position = [0.05 0.18 0.4 0.36];


%% Relative Stress Index (C)
load ./data/RF_v6_RSI_C_stats;

r = NaN(size(EcoL2));

% Add L1 stats
n = size(L1_stats, 1);
for j=1:n
    el = L1_stats(j, 1);
    rel = L1_stats(j, 3);

    if ~isnan(rel)
        r(EcoL1==el) = rel;
    end
end

% Add L2 stats
n = size(L2_stats, 1);
for j=1:n
    el = L2_stats(j, 1);
    rel = L2_stats(j, 3);

    if ~isnan(rel)
        r(EcoL2==el) = rel;
    end
end

% Add L3 stats
n = size(L3_stats, 1);
for j=1:n
    el = L3_stats(j, 1);
    rel = L3_stats(j, 3);

    if ~isnan(rel)
        r(EcoL3==el) = rel;
    end
end

subplot(2,2,4)
ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
    'on','PLineLocation',10,'MLineLocation',10,'MeridianLabel','off',...
    'ParallelLabel','off','GLineWidth',0.3,'Frame','off','FFaceColor',...
    'none', 'FontName', 'Times New Roman',...
    'FLineWidth',1, 'FontSize',7, 'MLabelParallel',25.01,...
    'GColor',[0.8 0.8 0.8], 'FontColor',[0.4 0.4 0.4]);
surfm(PRISMlat, PRISMlon, r.*CONUS);
caxis([0 1]);
colormap(cbrew);
geoshow(states,'FaceColor','none','LineWidth',0.3, 'EdgeColor',[0.4 0.4 0.4])
axis off;
axis image;
text(-0.42, 0.85, 'D', 'FontName','Times New Roman', 'FontSize',12, 'FontWeight','bold');
ax.Position = [0.5703 0.18 0.4 0.36];

h = colorbar('southoutside');
h.Position = [0.11 0.12 0.8 0.03];
h.FontName = 'Times New Roman';
h.FontSize = 10;
h.Ticks = 0:0.125:1;
h.TickLabels = {'0','','0.25','','0.5','','0.75','','1'};
h.TickLength = 0.025;
h.Label.String = 'R^{2}';

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/rf-model-skill-maps.tif')
close all;

