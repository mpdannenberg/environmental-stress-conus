% Make map of number of years in tree-ring stress index database

syear = 1971;
load ./data/TREESI_wEnvFactors;

n = [TREESI.END] - syear + 1;
lat = [TREESI.LAT];
lon = [TREESI.LON];

states = shaperead('usastatehi','UseGeoCoords', true);
latlim = [24 49];
lonlim = [-125 -67];

clr = wesanderson('fantasticfox1');
clr2 = make_cmap([1 1 1; sqrt(clr(4,:)); clr(4,:)], 8);

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 7 4];

ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',[28,34,40,46],'MLineLocation',10,'MeridianLabel','on',...
        'ParallelLabel','on','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Times New Roman', 'GColor',[0.8 0.8 0.8],...
        'FLineWidth',1, 'FontColor',[0.4 0.4 0.4], 'MLabelLocation',20,...
        'MLabelParallel',24.01);
scatterm(lat, lon, 25, n, 'filled', 'MarkerEdgeColor','k');
geoshow(states,'FaceColor','none','LineWidth',0.4)
caxis([10 45])
colormap(clr2(2:end, :))
axis off;
axis image;

cb = colorbar('eastoutside');
cb.TickLength = 0.08;
ylabel(cb, {'Number of years per site during', 'the study period (1971-present)'}, 'FontSize',11);

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/supplemental-sample-depth-per-site.tif')
close all;

