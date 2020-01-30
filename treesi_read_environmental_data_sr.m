% Add climate, topography, and soils to ITRDB relative stress structure

tic

load ITRDB;
TREESI = ITRDB_TW; clearvars -except TREESI;

TREESI(isnan([TREESI.LAT]) | isnan([TREESI.LON])) = []; % Get rid of sites that are missing coordinates
TREESI = TREESI([TREESI.LAT] > 28 & [TREESI.LAT] <= 49 & [TREESI.LON] > -125  & [TREESI.LON] < -66);

site = {TREESI.SITE};
TREESI = TREESI(~contains(site, 'cana')); % Remove excess Canadian sites
site = {TREESI.SITE};
TREESI = TREESI(~contains(site, 'mexi')); % Remove excess Mexico sites

% Load coordinates and convert to lat/lon array for calculating shortest
% distance
n = length(TREESI);
cd PRISM;
load prism_latlon;
nx = length(PRISMlon);
ny = length(PRISMlat);
prismLatLon = [reshape(repmat(PRISMlat', 1, nx), [], 1) reshape(repmat(PRISMlon, ny, 1), [], 1)];
syear = 1970;
eyear = 2013;

%% Read in EcoRegions
% Level 1
EcoL1 = double(geotiffread('./data/us_EcoL1_4km.tif'));
EcoL1(EcoL1==255) = NaN;

% Level 2
EcoL2 = double(geotiffread('./data/us_EcoL2_4km.tif'));
EcoL2(EcoL2==255) = NaN;

% Level 3
EcoL3 = double(geotiffread('./data/us_EcoL3_4km.tif'));
EcoL3(EcoL3==65536) = NaN;

%% Read in topographic data
% Elevation
ELEV = double(geotiffread('./data/us_dem_4km.tif'));
ELEV(ELEV<-1000) = NaN;

% Slope
SLOPE = double(geotiffread('./data/us_slope_4km.tif')) * 0.01;
SLOPE(SLOPE<-1000) = NaN;

% UAA (log-transformed)
UAA = double(geotiffread('./data/us_uaa_4km.tif'));
UAA(UAA<-1000) = NaN;
UAA = log(UAA);

% TWI 
TWI = double(geotiffread('./data/us_twi_4km.tif')) * 0.01;
TWI(TWI<-1000) = NaN;

%% Read in soil data (averaged for top 2.3 m of soil)
% Clay content (% of weight)
CLAY = double(geotiffread('./data/us_CLAY_4km.tif'));
CLAY(CLAY<-1000) = NaN;

% Organic carbon (% of weight)
OC = double(geotiffread('./data/us_OC_4km.tif'));
OC(OC<-1000) = NaN;

% pH (CaCl2 method)
PHCA = double(geotiffread('./data/us_PHCA_4km.tif'));
PHCA(PHCA<-1000) = NaN;

% Sand content (% of weight)
SAND = double(geotiffread('./data/us_SAND_4km.tif'));
SAND(SAND<-1000) = NaN;

% Silt content (% of weight)
SILT = double(geotiffread('./data/us_SILT_4km.tif'));
SILT(SILT<-1000) = NaN;

% Total N (% of weight)
TN = double(geotiffread('./data/us_TN_4km.tif'));
TN(TN<-1000) = NaN;

%% Loop through sites and extract environmental information
cd ../PRISM;
f1 = matfile('PRISM_PPT');
f2 = matfile('PRISM_TMAX');
f3 = matfile('PRISM_TMIN');
f4 = matfile('PRISM_PET.mat');
f5 = matfile('PRISM_VPD.mat');

for i = 1:n
    
    xy = [TREESI(i).LAT TREESI(i).LON];
    DistDeg = distance(xy(1), xy(2), prismLatLon(:,1), prismLatLon(:,2));
    DistKM = distdim(DistDeg, 'deg', 'km');
    
    xy = prismLatLon(DistKM == min(DistKM), :);
    xind = find(PRISMlon == xy(2));
    yind = find(PRISMlat == xy(1));
    
    % Ecoregions
    TREESI(i).EcoL1 = EcoL1(yind, xind);
    TREESI(i).EcoL2 = EcoL2(yind, xind);
    TREESI(i).EcoL3 = EcoL3(yind, xind);
    
    % Topography
    TREESI(i).ELEV = ELEV(yind, xind);
    TREESI(i).SLOPE = SLOPE(yind, xind);
    TREESI(i).UAA = UAA(yind, xind);
    TREESI(i).TWI = TWI(yind, xind);
    
    % Soils
    TREESI(i).GSDE.CLAY = CLAY(yind, xind);
    TREESI(i).GSDE.OC = OC(yind, xind);
    TREESI(i).GSDE.PHCA = PHCA(yind, xind);
    TREESI(i).GSDE.SAND = SAND(yind, xind);
    TREESI(i).GSDE.SILT = SILT(yind, xind);
    TREESI(i).GSDE.TN = TN(yind, xind);
    
    % Climate
    TREESI(i).PRISM.YEAR = syear:eyear;
    TREESI(i).PRISM.PPT = squeeze(f1.PPT(yind, xind, :, :));
    TREESI(i).PRISM.TMAX = squeeze(f2.TMAX(yind, xind, :, :));
    TREESI(i).PRISM.TMIN = squeeze(f3.TMIN(yind, xind, :, :));
    TREESI(i).PRISM.PET = squeeze(f4.PET(yind, xind, :, :));
    TREESI(i).PRISM.VPD = squeeze(f5.VPD(yind, xind, :, :));
    
    clear xy DistDeg DistKM xind yind;
    
end
clear f1 f2 f3 f4 f5 prismLatLon nx ny i syear eyear n;

cd ..

TREESI(isnan([TREESI.EcoL1])) = [];
save('TREESI_wEnvFactors_RSI','TREESI', '-v7.3');

toc

