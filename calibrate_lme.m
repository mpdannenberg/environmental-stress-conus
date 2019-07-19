%% Mixed effects models of tree-ring stress index
% M Dannenberg, 15 Apr 2016

%% Parameters
n_min = 75;
syear = 1971;

% Number of vars for each environmental stress type
nTopo = 4;
nSoil = 7; % Subject to change
nTemp = 8;
nWater = 8;

% Names of variables
TopoLabs = {'ELEV','SLOPE','UAA','TWI'};
SoilLabs = {'CLAY','OC','PHCA','SAND','SILT','TN','VMC2'};
TempLabs = {'TMIN_SON','TMIN_DJF','TMIN_MAM','TMIN_JJA',...
	'TMAX_SON','TMAX_DJF','TMAX_MAM','TMAX_JJA'};
WaterLabs = {'WB_SON','WB_DJF','WB_MAM','WB_JJA',...
	'VPD_SON','VPD_DJF','VPD_MAM','VPD_JJA'};
AllLabs = ['Stress','Year',TopoLabs SoilLabs TempLabs WaterLabs,'Site','Species'];

p = length(AllLabs);

% month of previous year to start with
smonth = 9; % September?

%% Load tree-ring stress estimates and put in table
load ./data/TREESI_wEnvFactors;

n = sum([TREESI.END] - syear + 1);

T = table('size',[n 2+nTopo+nSoil+nTemp+nWater+2], 'VariableTypes',...
    {'double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double',...
    'double','double','double','double','double','double','double','double','string','string'},...
    'VariableNames',AllLabs);
idx0 = 1;

for i = 1:length(TREESI)
    
    yr = TREESI(i).YEAR;
    ind = yr >= syear;
    idx = idx0:(idx0+sum(ind)-1);
    idx0 = idx0+sum(ind);
    
    T.Stress(idx) = TREESI(i).TREESI(ind);
    T.Year(idx) = yr(ind);
    T.ELEV(idx) = TREESI(i).ELEV;
    T.SLOPE(idx) = TREESI(i).SLOPE;
    T.UAA(idx) = TREESI(i).UAA;
    T.TWI(idx) = TREESI(i).TWI;
    T.CLAY(idx) = TREESI(i).GSDE.CLAY;
    T.OC(idx) = TREESI(i).GSDE.OC;
    T.PHCA(idx) = TREESI(i).GSDE.PHCA;
    T.SAND(idx) = TREESI(i).GSDE.SAND;
    T.SILT(idx) = TREESI(i).GSDE.SILT;
    T.TN(idx) = TREESI(i).GSDE.TN;
    T.VMC2(idx) = TREESI(i).GSDE.VMC2; % volumetric water content at -33 kPa

    T.Site(idx) = TREESI(i).SITE;
    T.Species(idx) = TREESI(i).SPECIES;
    
    [~, ind] = intersect(TREESI(i).PRISM.YEAR, yr(ind));
    % Ugh... gotta do all the calculations here...
    
    % 3-month P-PET sums
    P = reshape((TREESI(i).PRISM.PPT - TREESI(i).PRISM.PET)', 1, []);
    windowSize = 3;
    b = ones(1,windowSize);
    a = 1;
    P = filter(b, a, P);
    P = reshape(P, 12, [])';
    T.WB_SON(idx) = P(ind-1, 11);
    T.WB_DJF(idx) = P(ind, 2);
    T.WB_MAM(idx) = P(ind, 5);
    T.WB_JJA(idx) = P(ind, 8);
    
    % 3-month Tmin, Tmax, and VPD means
    Tmin = reshape(TREESI(i).PRISM.TMIN', 1, []);
    Tmax = reshape(TREESI(i).PRISM.TMAX', 1, []);
    VPD = reshape(TREESI(i).PRISM.VPD', 1, []);
    windowSize = 3;
    b = ones(1,windowSize)/windowSize;
    a = 1;
    Tmin = reshape(filter(b, a, Tmin), 12, [])';
    Tmax = reshape(filter(b, a, Tmax), 12, [])';
    VPD = reshape(filter(b, a, VPD), 12, [])';
    T.TMIN_SON(idx) = Tmin(ind-1, 11);
    T.TMIN_DJF(idx) = Tmin(ind, 2);
    T.TMIN_MAM(idx) = Tmin(ind, 5);
    T.TMIN_JJA(idx) = Tmin(ind, 8);
    T.TMAX_SON(idx) = Tmax(ind-1, 11);
    T.TMAX_DJF(idx) = Tmax(ind, 2);
    T.TMAX_MAM(idx) = Tmax(ind, 5);
    T.TMAX_JJA(idx) = Tmax(ind, 8);
    T.VPD_SON(idx) = VPD(ind-1, 11);
    T.VPD_DJF(idx) = VPD(ind, 2);
    T.VPD_MAM(idx) = VPD(ind, 5);
    T.VPD_JJA(idx) = VPD(ind, 8);
    
    
end
writetable(T, 'StressIndexCovariateTable.csv');

clearvars -except T;

%% Calibrate LME
FE_TMIN = 'TMIN_SON + TMIN_DJF + TMIN_MAM + TMIN_JJA';
FE_TMAX = 'TMAX_SON + TMAX_DJF + TMAX_MAM + TMAX_JJA';
FE_WB = 'WB_SON + WB_DJF + WB_MAM + WB_JJA'; 
FE_VPD = 'VPD_SON + VPD_DJF + VPD_MAM + VPD_JJA';
FE_Site = '(ELEV + SLOPE + UAA + TWI + SAND + SILT + CLAY + OC + TN + PHCA + VMC2)*';

% Initial full model - selection by AIC
mdl = ['Stress ~ ', FE_Site,'(',FE_TMIN,'+',FE_TMAX,'+',FE_WB,'+',FE_VPD,')'];
lme0 = fitlme(T, [mdl,'+ (1 | Year) + (1 | Species) + (1 | Species:Site)']);

aic0 = 0;
aic1 = lme0.ModelCriterion.AIC;

while aic1 < aic0
    
    aic0 = aic1;
    lme0 = lme1;
    
    p = lme0.Coefficients.pValue(2:end);
    var = lme0.Coefficients.Name(2:end);
    
    var2remove = var{p == max(p)};
    mdl = [mdl, ' - ',var2remove];
    
    lme1 = fitlme(T, [mdl,'+ (1 | Year) + (1 | Species) + (1 | Species:Site)']);
    
    aic1 = lme1.ModelCriterion.AIC;
    
end

% Initial full model - selection by BIC
mdl = ['Stress ~ ', FE_Site,'(',FE_TMIN,'+',FE_TMAX,'+',FE_WB,'+',FE_VPD,')'];
lme0 = fitlme(T, [mdl,'+ (1 | Year) + (1 | Species) + (1 | Species:Site)']);

bic0 = 0;
bic1 = lme0.ModelCriterion.BIC;

while bic1 < bic0
    
    bic0 = bic1;
    lme0 = lme1;
    
    p = lme0.Coefficients.pValue(2:end);
    var = lme0.Coefficients.Name(2:end);
    
    var2remove = var{p == max(p)};
    mdl = [mdl, ' - ',var2remove];
    
    lme1 = fitlme(T, [mdl,'+ (1 | Year) + (1 | Species) + (1 | Species:Site)']);
    
    bic1 = lme1.ModelCriterion.BIC;
    
end



