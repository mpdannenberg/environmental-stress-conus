% Map L1 ecoregions w/ tree-ring sites and AmeriFlux sites
% Changes needed:
    % 

load AmeriFlux;
load TREESI_wEnvFactors;

TLAT = [TREESI.LAT];
TLON = [TREESI.LON];
TSTART = [TREESI.START];
TEND = [TREESI.END];

AmeriFlux = AmeriFlux(~cellfun(@isempty,{AmeriFlux.LAT}));
SITES = {AmeriFlux.SITE};

StressFluxComp = NaN(1000, 12); % Column order: Stress, GPP, GPP anomaly, distance, elev difference
StressFluxComp_info = cell(1000, 3); % Column order: site, species
sfidx = 1;

for i = 1:length(AmeriFlux)
    FluxDates = AmeriFlux(i).YEAR;
    FluxLat = AmeriFlux(i).LAT;
    FluxLon = AmeriFlux(i).LON;
    gpp = AmeriFlux(i).ANNUAL.GPP;
    lue = AmeriFlux(i).ANNUAL.LUE;
    DistDeg = distance(FluxLat, FluxLon, TLAT, TLON);
    DistKM = distdim(DistDeg, 'deg', 'km');
    
    NumWithin100km(i) = sum(DistKM<=100);
    inds = find(DistKM<=100);
    dist = DistKM(inds);
    if sum(DistKM<=100)>=1
        for j = 1:length(inds)
            idx = inds(j);
            TreeDates = TREESI(idx).YEAR;
            
            temp(j) = sum(ismember(TreeDates,FluxDates(~isnan(lue))));
            
            if max(temp(j)) > 0
                
                
                tidx = ismember(TreeDates,FluxDates);
                fidx = ismember(FluxDates,TreeDates);
                stress = TREESI(idx).TREESI(tidx);
                n = length(stress);
                StressFluxComp(sfidx:(sfidx+n-1), 1) = stress;
                StressFluxComp(sfidx:(sfidx+n-1), 2) = gpp(fidx);
                StressFluxComp(sfidx:(sfidx+n-1), 3) = gpp(fidx) - nanmean(gpp);
                StressFluxComp(sfidx:(sfidx+n-1), 4) = lue(fidx);
                StressFluxComp(sfidx:(sfidx+n-1), 5) = dist(j);
                StressFluxComp(sfidx:(sfidx+n-1), 6) = abs( AmeriFlux(i).ELEV - TREESI(idx).ELEV );
                StressFluxComp(sfidx:(sfidx+n-1), 7) = AmeriFlux(i).LAT;
                StressFluxComp(sfidx:(sfidx+n-1), 8) = AmeriFlux(i).LON;
                StressFluxComp(sfidx:(sfidx+n-1), 9) = TREESI(idx).LAT;
                StressFluxComp(sfidx:(sfidx+n-1), 10) = TREESI(idx).LON;
                StressFluxComp(sfidx:(sfidx+n-1), 11) = AmeriFlux(i).ANNUAL.MOD17.fE(fidx);
                StressFluxComp(sfidx:(sfidx+n-1), 12) = AmeriFlux(i).ANNUAL.ECLUE.fE(fidx);
                for ThisCell = 1:n
                    StressFluxComp_info{(sfidx+ThisCell-1), 1} = TREESI(idx).SITE;
                    StressFluxComp_info{(sfidx+ThisCell-1), 2} = TREESI(idx).SPECIES;
                    StressFluxComp_info{(sfidx+ThisCell-1), 3} = AmeriFlux(i).IGBP;
                end
                sfidx = sfidx + n;
                    

            end
        end
        MaxOverlap(i) = max(max(temp));
        
        
    else
        MaxOverlap(i) = 0;
    end
    
    
end

StressFluxComp_info = StressFluxComp_info(~any(isnan(StressFluxComp(:,1:3)),2),:);
StressFluxComp = StressFluxComp(~any(isnan(StressFluxComp(:,1:3)),2),:);


%load TREESI_wEnvFactors3;
load conus_mask;
cd EcoRegions;
ecol1 = double(geotiffread('us_EcoL1_4km.tif'));
ecol1(ecol1==255 | ecol1==0 | ecol1==150) = NaN;
ecol1 = ecol1.*CONUS;
cd ..

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

% cbrew = [31,120,180
%     178,223,138
%     51,160,44
%     166,206,227
%     251,154,153
%     253,191,111
%     255,127,0
%     202,178,214
%     106,61,154]/255;

% cbrew = [190,186,218
% 255,255,179
% 141,211,199
% 251,128,114
% 128,177,211
% 253,180,98
% 179,222,105
% 252,205,229
% 217,217,217]/255;

cbrew = [251,180,174
179,205,227
204,235,197
222,203,228
254,217,166
255,255,204
229,216,189
253,218,236
232,232,232]/255;

% n = length(AmeriFlux);
% nyrs = NaN(1,n);
% for i = 1:n
%    
%     nyrs(i) = sum(isfinite(AmeriFlux(i).ANNUAL.GPP));
%     
% end
% AmeriFlux = AmeriFlux(nyrs > 0);

% alat = [AmeriFlux.LAT];
% alon = [AmeriFlux.LON];
alat = StressFluxComp(:,7);
alon = StressFluxComp(:,8);

tlat = [TREESI.LAT];
tlon = [TREESI.LON];

load TREESI_wEnvFactors3;
tlat2 = [TREESI.LAT];
tlon2 = [TREESI.LON];




% Map time
cd PRISM

load prism_latlon;
cd ..;

states = shaperead('usastatehi','UseGeoCoords', true);
% worldland= shaperead('landareas', 'UseGeoCoords', true);
latlim = [24 50];
lonlim = [-125 -67];
tr_lat = [46.6991005000000;48.7808220000000;48.1821765000000;48.5064070000000;48.5867320000000;48.5967320000000;47.3394060000000];
tr_lon = [-120.932076500000;-119.272581500000;-120.253510000000;-118.754520000000;-119.138224991000;-119.696205000000;-120.554938000000];

h=figure('Color','w');
h.Units = 'inches';
h.Position = [1 4 7.3 3.5];
% set(h, 'Position',[2 42 638 594]);

subplot(1,3,[1 2])
ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',10,'MLineLocation',10,'MeridianLabel','on',...
        'ParallelLabel','on','GLineWidth',0.3,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Times New Roman', 'GColor',[0.8 0.8 0.8],...
        'FLineWidth',1, 'FontColor',[0.4 0.4 0.4], 'MLabelLocation',20,...
        'MLabelParallel',24.01);
surfm(PRISMlat, PRISMlon, ecol1_reclass);
geoshow(states,'FaceColor','none','LineWidth',0.4)
% p3 = plotm(tlat2, tlon2, 'k^','MarkerFaceColor','w','MarkerSize',2);
p3 = plotm(tlat2, tlon2, 'k.','MarkerSize',6);
p1 = plotm(tlat, tlon, 'k^','MarkerFaceColor','w','MarkerSize',4, 'LineWidth',0.3);
plotm(tr_lat, tr_lon, 'k^','MarkerFaceColor','w','MarkerSize',4, 'LineWidth',0.3);
p2 = plotm(alat, alon, 'ko','MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',6, 'LineWidth',1.2);
caxis([0.5 9.5])
colormap(cbrew)
axis off;
axis image;
set(ax, 'Position',[0.02 0.1 0.82 0.87]);

% annotation('rectangle',[0.7 0 0.3 1],'Color','w', 'FaceColor','w')

h = colorbar('eastoutside');
pos = get(h, 'Position');
set(h, 'Position',[0.73 0.15 0.0381 0.5876],'FontName','Times New Roman',...
    'XTick',1:1:9, 'XTickLabel',ecolabs, 'FontSize',8);

h = legend([p3,p1,p2],'Tree-ring site ({\itS_{r}} only)','Tree-ring site ({\itS^{*}} and {\itS_{r}})','Flux Tower');
set(h, 'Position',[0.73 0.77 0.1982 0.18], 'FontName','Times New Roman', 'FontSize',8);
legend('boxoff');

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r600','EcoL1_ITRDB_AmeriFlux_map4.tif')



% Note: could look into new flux dataset
