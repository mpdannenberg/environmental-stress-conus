% Map L1 ecoregions w/ tree-ring sites and AmeriFlux sites

load ./data/TREESI_wEnvFactors;

TLAT = [TREESI.LAT];
TLON = [TREESI.LON];
TSTART = [TREESI.START];
TEND = [TREESI.END];

load ./data/conus_mask;
ecol1 = double(geotiffread('./data/us_EcoL1_4km.tif'));
ecol1(ecol1==255 | ecol1==0 | ecol1==150) = NaN;
ecol1 = ecol1.*CONUS;

clear CONUS;

ecos = unique( ecol1(isfinite(ecol1))' );

ecol1_reclass = NaN(size(ecol1));
for i=1:length(ecos)
    
    ecol1_reclass(ecol1 == ecos(i)) = i;
    
end

ecolabs = {'Northern Forests','Northwestern Forested Mountains',...
    'Marine West Coast Forest','Eastern Temperate Forests','Great Plains',...
    'North American Deserts','Mediterranean California',...
    'Southern Semi-Arid Highlands','Temperate Sierras'};
ecolabs2 = {'5.0','6.0','7.0','8.0','9.0','10.0','11.0','12.0','13.0'};

clr = wesanderson('fantasticfox1');
cbrew = [clr(4,:)
    sqrt(clr(4,:))
    clr(3,:).^2
    clr(3,:)
    clr(1,:).^2
    clr(1,:)
    clr(2,:)
    sqrt(clr(2,:))
    clr(5,:)];

tlat = [TREESI.LAT];
tlon = [TREESI.LON];

load ./data/TREESI_wEnvFactors_RSI;
tlat2 = [TREESI.LAT];
tlon2 = [TREESI.LON];

% Time to make a map...
load ./data/prism_latlon;

states = shaperead('usastatehi','UseGeoCoords', true);
latlim = [24 49];
lonlim = [-125 -67];
tr_lat = [46.6991005000000;48.7808220000000;48.1821765000000;48.5064070000000;48.5867320000000;48.5967320000000;47.3394060000000];
tr_lon = [-120.932076500000;-119.272581500000;-120.253510000000;-118.754520000000;-119.138224991000;-119.696205000000;-120.554938000000];

h=figure('Color','w');
h.Units = 'inches';
h.Position = [1 4 6.5 3];

subplot(1,3,[1 2])
ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',[28,34,40,46],'MLineLocation',10,'MeridianLabel','on',...
        'ParallelLabel','on','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Times New Roman', 'GColor',[0.8 0.8 0.8],...
        'FLineWidth',1, 'FontColor',[0.4 0.4 0.4], 'MLabelLocation',20,...
        'MLabelParallel',24.01);
surfm(PRISMlat, PRISMlon, ecol1_reclass);
geoshow(states,'FaceColor','none','LineWidth',0.4)
p3 = plotm(tlat2, tlon2, 'k.','MarkerSize',6);
p1 = plotm(tlat, tlon, 'k^','MarkerFaceColor','w','MarkerSize',4, 'LineWidth',0.3);
plotm(tr_lat, tr_lon, 'k^','MarkerFaceColor','w','MarkerSize',4, 'LineWidth',0.3);
caxis([0.5 9.5])
colormap(cbrew)
axis off;
axis image;
set(ax, 'Position',[0.03 0.08 0.82 0.87]);

h = colorbar('eastoutside');
pos = get(h, 'Position');
set(h, 'Position',[0.71 0.1 0.0381 0.5876],'FontName','Times New Roman',...
    'XTick',1:1:9, 'XTickLabel',ecolabs, 'FontSize',8);

h = legend([p3,p1],'{\itS_{r}} only','{\itS^{*}} and {\itS_{r}}');
set(h, 'Position',[0.7 0.75 0.1982 0.18], 'FontName','Times New Roman', 'FontSize',8);
legend('boxoff');

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/fig1-site-map.tif')
close all;


