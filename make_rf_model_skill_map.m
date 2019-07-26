%%
%% Map model strength for TSC and C models for both SIa and SIr
%%


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


% Order of *_stats columns: ecoregion, rmse, r2, Spearman's rho, n, nsites

% cbrew = [247,252,240
% 224,243,219
% 204,235,197
% 168,221,181
% 123,204,196
% 78,179,211
% 43,140,190
% 8,88,158]/255;

cbrew = [255,247,251
236,231,242
208,209,230
166,189,219
116,169,207
54,144,192
5,112,176
3,78,123]/255;

h=figure('Color','w');
h.Units = 'inches';
h.Position = [2 4 6.65 4];

%% Absolute Stress Index (TSC)
cd RF_v4_TSC;
load RF_v4_TSC_stats;
cd ..;

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
% title('SI_{a,TSC}', 'FontName','Times New Roman');
text(-0.42, 0.85, '(a)', 'FontName','Times New Roman', 'FontSize',12, 'FontWeight','bold');
text(-0.05, 0.96, 'TSC model', 'FontName','Times New Roman', 'FontSize',14,...
    'FontWeight','bold','HorizontalAlignment','center');
ax.Position = [0.05 0.55 0.4 0.36];

%% Absolute Stress Index (C)
cd RF_v4_C;
load RF_v4_C_stats;
cd ..;

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
% title('SI_{a,C}', 'FontName','Times New Roman');
text(-0.42, 0.85, '(b)', 'FontName','Times New Roman', 'FontSize',12, 'FontWeight','bold');
text(-0.05, 0.96, 'C model', 'FontName','Times New Roman', 'FontSize',14,...
    'FontWeight','bold','HorizontalAlignment','center');
ax.Position = [0.5703 0.55 0.4 0.36];


%% Relative Stress Index (TSC)
cd RF_v6_RSI_TSC;
load RF_v6_RSI_TSC_stats;
cd ..;

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
% title('SI_{r,TSC}', 'FontName','Times New Roman');
text(-0.42, 0.85, '(c)', 'FontName','Times New Roman', 'FontSize',12, 'FontWeight','bold');
ax.Position = [0.05 0.05 0.4 0.36];


%% Relative Stress Index (C)
cd RF_v6_RSI_C;
load RF_v6_RSI_C_stats;
cd ..;

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
% title('SI_{r,C}', 'FontName','Times New Roman');
text(-0.42, 0.85, '(d)', 'FontName','Times New Roman', 'FontSize',12, 'FontWeight','bold');
ax.Position = [0.5703 0.05 0.4 0.36];

h = colorbar('southoutside');
h.Position = [0.44 0.12 0.15 0.03];
h.FontName = 'Times New Roman';
h.FontSize = 10;
h.Ticks = 0:0.25:1;
h.Label.String = 'r^{2}';

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r600','RF_TSCvsC_ModelStrength.tif')
close all;

