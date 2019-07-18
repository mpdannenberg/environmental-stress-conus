rs_yr = 1982:2015;

% Flux data
[num, txt] = xlsread('us-locs-20171009140829.xlsx',1, 'A1:E64');

flx_lat = num(:,1);
flx_lon = num(:,2);
clear num txt;

[num, txt] = xlsread('us-availability-20171009140829.xlsx',1, 'B2:Y65');
flx_t = sum(num>0);
flx_x = sum(num>0, 2);
flx_yr = 1991:2014;
clear num txt;


% Tree rings
load ITRDB;
ITRDB = ITRDB_TW;
clear ITRDB_EW ITRDB_LWadj ITRDB_TW;
ITRDB = ITRDB([ITRDB.LAT] > 28 & [ITRDB.LAT] <= 49 & [ITRDB.LON] > -125  & [ITRDB.LON] < -66);
ITRDB(373:395) = []; % get rid of excess Canada
ITRDB(656:673) = []; % get rid of excess Mexico
ITRDB(1277) = []; % get rid of silly one

tr_lat = [ITRDB.LAT];
tr_lon = [ITRDB.LON];

num = zeros(length(ITRDB), length(rs_yr));
for i = 1:length(ITRDB)
    
    year = ITRDB(i).YEAR;
    
    for j = 1:length(rs_yr)
        
        idx = year == rs_yr(j);
        num(i, j) = sum(idx);
        
    end
    
end

tr_t = sum(num>0);
tr_x = sum(num>0, 2);
clear num i j idx year;

% RS overlap
rs_t = zeros(2, length(rs_yr));
rs_t(1, :) = tr_t;
[~, idx, ~] = intersect(rs_yr, flx_yr);
rs_t(2, idx) = flx_t;
clear idx;

% % Plot
% h = figure('Color','w');
% h.Units = 'inches';
% h.Position = [1 1 6 4];
% 
% plot(rs_yr, tr_t, '-', 'LineWidth',2, 'Color',[51,160,44]/255)
% set(gca, 'XLim',[1982 2015])
% hold on;
% plot(flx_yr, flx_t, '--', 'LineWidth',2, 'Color',[31,120,180]/255)
% close all;

% Alt plot
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6 4];

patch([rs_yr 2015 1982], [tr_t 0 0], [51,160,44]/255)
set(gca, 'XLim',[1982 2015])
set(gca, 'FontSize',12);alpha(0.5);
hold on;
patch([flx_yr 2015 min(flx_yr)], [flx_t 0 0], [31,120,180]/255)
alpha(0.5)
grid on;
box on;
[BL,BLicons] = legend('ITRDB','FLUXNET 2015');
legend('boxoff');
PatchInLegend = findobj(BLicons, 'type', 'patch');
set(PatchInLegend, 'facea', 0.5)

xlabel('Year', 'FontSize',16);
ylabel('Number of sites', 'FontSize',16);

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r600','FlxTR_siteyrs.tif')
close all;

% Map
states = shaperead('usastatehi','UseGeoCoords', true);
latlim = [25.1 49.5];
lonlim = [-125 -66];

h=figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 7 4];

ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',5,'MLineLocation',10,'MeridianLabel','on',...
        'ParallelLabel','on','GLineWidth',0.5,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',25.11);
geoshow(states,'FaceColor',[0.95 0.95 0.95],'LineWidth',0.4)
axis off;
axis image;
p1 = scatterm(tr_lat(tr_x>0), tr_lon(tr_x>0), 8,'^', 'MarkerFaceColor',[178,223,138]/255, 'MarkerEdgeColor',[0.3 0.3 0.3]);
p2 = scatterm(flx_lat, flx_lon, 20, 'o','MarkerFaceColor',[166,206,227]/255,'MarkerEdgeColor',[0.3 0.3 0.3], 'LineWidth',1.2);

h = legend([p1,p2],'ITRDB','FLUXNET 2015');
set(h, 'Position',[0.1 0.8 0.8 0.1302], 'FontName','Helvetica');
legend('boxoff');

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r600','Flx_sites.tif')
close all;



%% Now replace with TREESI

load TREESI_wEnvFactors.mat;

tr_lat = [TREESI.LAT];
tr_lon = [TREESI.LON];

num = zeros(length(TREESI), length(rs_yr));
for i = 1:length(TREESI)
    
    year = TREESI(i).YEAR;
    
    for j = 1:length(rs_yr)
        
        idx = year == rs_yr(j);
        num(i, j) = sum(idx);
        
    end
    
end

tr_t = sum(num>0);
tr_x = sum(num>0, 2);
clear num i j idx year;

% Alt plot
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6 4];

patch([rs_yr 2015 1982], [tr_t 0 0], [51,160,44]/255)
set(gca, 'XLim',[1982 2015])
set(gca, 'FontSize',12);alpha(0.5);
hold on;
patch([flx_yr 2015 min(flx_yr)], [flx_t 0 0], [31,120,180]/255)
alpha(0.5)
grid on;
box on;
[BL,BLicons] = legend('ITRDB Stress Index','FLUXNET 2015');
legend('boxoff');
PatchInLegend = findobj(BLicons, 'type', 'patch');
set(PatchInLegend, 'facea', 0.5)

xlabel('Year', 'FontSize',16);
ylabel('Number of sites', 'FontSize',16);

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r600','FlxTREESI_siteyrs.tif')
close all;

% Map

% load('C:\Users\Matt\Documents\Publications\Dannenberg_TR_NDVI\conus_mask.mat')
% load('C:\Users\Matt\Documents\Publications\Dannenberg_TR_NDVI\prism_latlon.mat')
% cd('C:\Users\Matt\Desktop\Dissertation\Chpt3\EcoRegions');
% ecol1 = double(geotiffread('us_EcoL1_4km.tif'));
% ecol1(ecol1==255 | ecol1==0 | ecol1==150) = NaN;
% ecol1 = ecol1.*CONUS;
% cd('C:\Users\Matt\Documents\Publications\Dannenberg_et_al_GlobalEcolBiogeog');
% 
% clear CONUS;
% 
% ecos = unique( ecol1(isfinite(ecol1))' );
% ecol1_reclass = NaN(size(ecol1));
% for i=1:length(ecos)
%     
%     ecol1_reclass(ecol1 == ecos(i)) = i;
%     
% end
% 
% ecolabs = {'Northern Forests','Northwestern Forested Mountains',...
%     'Marine West Coast Forest','Eastern Temperate Forests','Great Plains',...
%     'North American Deserts','Mediterranean California',...
%     'Southern Semi-Arid Highlands','Temperate Sierras'};
% ecolabs2 = {'5.0','6.0','7.0','8.0','9.0','10.0','11.0','12.0','13.0'};
% 
% cbrew = [31,120,180
%     178,223,138
%     51,160,44
%     166,206,227
%     251,154,153
%     253,191,111
%     255,127,0
%     202,178,214
%     106,61,154]/255;

states = shaperead('usastatehi','UseGeoCoords', true);
latlim = [25.1 49.5];
lonlim = [-125 -66];

h=figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 7 4];

ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',5,'MLineLocation',10,'MeridianLabel','on',...
        'ParallelLabel','on','GLineWidth',0.5,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',25.11);
% h=surfm(PRISMlat, PRISMlon, ecol1_reclass);
geoshow(states,'FaceColor',[0.95 0.95 0.95],'LineWidth',0.4)
axis off;
axis image;
p1 = scatterm(tr_lat, tr_lon, 12,'^', 'MarkerFaceColor',[178,223,138]/255, 'MarkerEdgeColor',[0.3 0.3 0.3]);
p2 = scatterm(flx_lat, flx_lon, 20, 'o','MarkerFaceColor',[166,206,227]/255,'MarkerEdgeColor',[0.3 0.3 0.3], 'LineWidth',1.2);
% caxis([0.5 9.5])
% colormap(cbrew)

h = legend([p1,p2],'ITRDB Stress Index','FLUXNET 2015');
set(h, 'Position',[0.1 0.82 0.8 0.1302], 'FontName','Helvetica');
legend('boxoff');

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r600','FlxTREESI_sites.tif')
close all;

