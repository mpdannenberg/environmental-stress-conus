% Read stress index files produced from ITRDB rwls
% M Dannenberg 8/9/14
% modified 9/12/15
% modified 9/20/16: read in for relative stress index

tic

% set directory to local ITRDB directory

[ninfo, tinfo, ~] = xlsread('ITRDB_metadata.xlsx','A:G');
SITE = tinfo(2:end, 1);
SPECIES = tinfo(2:end, 4);
LAT = ninfo(:, 4);
LON = ninfo(:, 5);
ELEV = ninfo(:, 6);
clear ninfo tinfo;

% TRW
fnames = glob('*w_crn.csv');

for i = 1:length(fnames)
    sname = regexp(fnames{i}, '\w_', 'split');
    TF = strcmpi(sname{1}, SITE);
    idx = find(TF);
    
    % Basic Info
    ITRDB_TW(i).SITE = SITE{idx};
    ITRDB_TW(i).SPECIES = SPECIES{idx};
    ITRDB_TW(i).LAT = LAT(idx);
    ITRDB_TW(i).LON = LON(idx);
    ITRDB_TW(i).ELEV = ELEV(idx);
    
    
    % Read CSV with RWI
    out = read_mixed_csv(fnames{i}, ',');
    std = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,2)));
    res = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,3)));
    samp = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,4)));
    year = str2double(strrep(out(2:end, 1), '"', ''));
    
    ITRDB_TW(i).STD = std;
    ITRDB_TW(i).RES = res;
    ITRDB_TW(i).SAMPLE_DEPTH = samp;
    ITRDB_TW(i).YEAR = year;
    ITRDB_TW(i).START = min(year);
    ITRDB_TW(i).END = max(year);
    
    clear out std res samp year;
    
    
    % Read CSV with stats
    out = read_mixed_csv([sname{1},'w_stats.csv'], ',');
    ITRDB_TW(i).STATS.StartYear = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,2)));
    ITRDB_TW(i).STATS.MidYear = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,3)));
    ITRDB_TW(i).STATS.EndYear = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,4)));
    ITRDB_TW(i).STATS.n_Cores = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,5)));
    ITRDB_TW(i).STATS.n_Trees = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,6)));
    ITRDB_TW(i).STATS.n = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,7)));
    ITRDB_TW(i).STATS.RBAR_total = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,11)));
    ITRDB_TW(i).STATS.RBAR_WithinTree = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,12)));
    ITRDB_TW(i).STATS.RBAR_BetweenTree = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,13)));
    ITRDB_TW(i).STATS.RBAR_Effective = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,15)));
    ITRDB_TW(i).STATS.EPS = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,16)));
    ITRDB_TW(i).STATS.SNR = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,17)));
    
    clear out;
    clear sname TF idx;
    
    if rem(i, 100)==0
        time = toc;
        disp([num2str(i),' sites out of ', num2str(length(fnames)), ' complete: ', num2str(time), ' seconds'])
        tic
    end
end


% EW
fnames = glob('*e_crn.csv');

for i = 1:length(fnames)
    sname = regexp(fnames{i}, '\e_', 'split');
    TF = strcmpi(sname{1}, SITE);
    idx = find(TF);
    if isempty(idx)
        TF = strcmpi([sname{1},'w'], SITE);
        idx = find(TF);
    end
    
    % Basic Info
    ITRDB_EW(i).SITE = SITE{idx};
    ITRDB_EW(i).SPECIES = SPECIES{idx};
    ITRDB_EW(i).LAT = LAT(idx);
    ITRDB_EW(i).LON = LON(idx);
    ITRDB_EW(i).ELEV = ELEV(idx);
    
    
    % Read CSV with RWI
    out = read_mixed_csv(fnames{i}, ',');
    std = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,2)));
    res = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,3)));
    samp = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,4)));
    year = str2double(strrep(out(2:end, 1), '"', ''));
    
    ITRDB_EW(i).STD = std;
    ITRDB_EW(i).RES = res;
    ITRDB_EW(i).SAMPLE_DEPTH = samp;
    ITRDB_EW(i).YEAR = year;
    ITRDB_EW(i).START = min(year);
    ITRDB_EW(i).END = max(year);
    
    clear out std res samp year;
    
    
    % Read CSV with stats
    out = read_mixed_csv([sname{1},'e_stats.csv'], ',');
    ITRDB_EW(i).STATS.StartYear = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,2)));
    ITRDB_EW(i).STATS.MidYear = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,3)));
    ITRDB_EW(i).STATS.EndYear = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,4)));
    ITRDB_EW(i).STATS.n_Cores = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,5)));
    ITRDB_EW(i).STATS.n_Trees = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,6)));
    ITRDB_EW(i).STATS.n = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,7)));
    ITRDB_EW(i).STATS.RBAR_total = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,11)));
    ITRDB_EW(i).STATS.RBAR_WithinTree = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,12)));
    ITRDB_EW(i).STATS.RBAR_BetweenTree = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,13)));
    ITRDB_EW(i).STATS.RBAR_Effective = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,15)));
    ITRDB_EW(i).STATS.EPS = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,16)));
    ITRDB_EW(i).STATS.SNR = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,17)));
    
    clear out;
    clear sname TF idx;
    
    if rem(i, 100)==0
        time = toc;
        disp([num2str(i),' sites out of ', num2str(length(fnames)), ' complete: ', num2str(time), ' seconds'])
        tic
    end
end


% LWadj
fnames = glob('*la_crn.csv');

for i = 1:length(fnames)
    sname = regexp(fnames{i}, '\la_', 'split');
    TF = strcmpi(sname{1}, SITE);
    idx = find(TF);
    if isempty(idx)
        TF = strcmpi([sname{1},'w'], SITE);
        idx = find(TF);
    end

    % Basic Info
    ITRDB_LWadj(i).SITE = SITE{idx};
    ITRDB_LWadj(i).SPECIES = SPECIES{idx};
    ITRDB_LWadj(i).LAT = LAT(idx);
    ITRDB_LWadj(i).LON = LON(idx);
    ITRDB_LWadj(i).ELEV = ELEV(idx);
    
    
    % Read CSV with RWI
    out = read_mixed_csv(fnames{i}, ',');
    std = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,2)));
    res = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,3)));
    samp = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,4)));
    year = str2double(strrep(out(2:end, 1), '"', ''));
    
    ITRDB_LWadj(i).STD = std;
    ITRDB_LWadj(i).RES = res;
    ITRDB_LWadj(i).SAMPLE_DEPTH = samp;
    ITRDB_LWadj(i).YEAR = year;
    ITRDB_LWadj(i).START = min(year);
    ITRDB_LWadj(i).END = max(year);
    
    clear out std res samp year;
    
    
    % Read CSV with stats
    try
        out = read_mixed_csv([sname{1},'la_stats.csv'], ',');
        ITRDB_LWadj(i).STATS.StartYear = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,2)));
        ITRDB_LWadj(i).STATS.MidYear = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,3)));
        ITRDB_LWadj(i).STATS.EndYear = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,4)));
        ITRDB_LWadj(i).STATS.n_Cores = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,5)));
        ITRDB_LWadj(i).STATS.n_Trees = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,6)));
        ITRDB_LWadj(i).STATS.n = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,7)));
        ITRDB_LWadj(i).STATS.RBAR_total = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,11)));
        ITRDB_LWadj(i).STATS.RBAR_WithinTree = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,12)));
        ITRDB_LWadj(i).STATS.RBAR_BetweenTree = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,13)));
        ITRDB_LWadj(i).STATS.RBAR_Effective = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,15)));
        ITRDB_LWadj(i).STATS.EPS = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,16)));
        ITRDB_LWadj(i).STATS.SNR = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,17)));
    catch
        ITRDB_LWadj(i).STATS = 'Unavailable: adjustment done on site-level EW/LW';
    end
    
    clear out;
    clear sname TF idx;
    
    if rem(i, 100)==0
        time = toc;
        disp([num2str(i),' sites out of ', num2str(length(fnames)), ' complete: ', num2str(time), ' seconds'])
        tic
    end
end


save('ITRDB.mat', 'ITRDB_TW', 'ITRDB_EW', 'ITRDB_LWadj');
clear ELEV fnames i LAT LON SITE SPECIES time;

toc

%% Map

% TW
LAT = [ITRDB_TW.LAT];
LON = [ITRDB_TW.LON];
states = shaperead('usastatehi','UseGeoCoords',true);
worldland= shaperead('landareas', 'UseGeoCoords', true);
latlim = [-60 80];
lonlim = [-180 180];
h = figure('Color','w');
set(h, 'Units','inches', 'Position',[1 1 6.5 3]);

axis off; % comment out if using usamap
ax = axesm('pcarree','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
    'on','PLineLocation',20,'MLineLocation',60,'MeridianLabel',...
    'on','ParallelLabel','on','GLineWidth',0.5,'Frame','off','FontName',...
    'Times New Roman', 'GColor',[0.5 0.5 0.5], 'GLineStyle',':',...
    'GLineWidth',0.3, 'FontSize',8);
box off; % comment out if using usamap
geoshow(worldland,'FaceColor',[0.9 0.9 0.9],'LineWidth',0.4,'EdgeColor',[0.4 0.4 0.4])
geoshow(states,'FaceColor',[0.9 0.9 0.9],'LineWidth',0.4,'EdgeColor',[0.4 0.4 0.4])
plotm(LAT, LON, 'k^', 'MarkerSize',2, 'MarkerFaceColor',[208,209,230]/255);
axis image;

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','ITRDB_TW.tif')
close all;

% EW
LAT = [ITRDB_EW.LAT];
LON = [ITRDB_EW.LON];
states = shaperead('usastatehi','UseGeoCoords',true);
worldland= shaperead('landareas', 'UseGeoCoords', true);
latlim = [-60 80];
lonlim = [-180 180];
h = figure('Color','w');
set(h, 'Units','inches', 'Position',[1 1 6.5 3]);

axis off; % comment out if using usamap
ax = axesm('pcarree','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
    'on','PLineLocation',20,'MLineLocation',60,'MeridianLabel',...
    'on','ParallelLabel','on','GLineWidth',0.5,'Frame','off','FontName',...
    'Times New Roman', 'GColor',[0.5 0.5 0.5], 'GLineStyle',':',...
    'GLineWidth',0.3, 'FontSize',8);
box off; % comment out if using usamap
geoshow(worldland,'FaceColor',[0.9 0.9 0.9],'LineWidth',0.4,'EdgeColor',[0.4 0.4 0.4])
geoshow(states,'FaceColor',[0.9 0.9 0.9],'LineWidth',0.4,'EdgeColor',[0.4 0.4 0.4])
plotm(LAT, LON, 'k^', 'MarkerSize',3, 'MarkerFaceColor',[214,96,77]/255);
LAT = [ITRDB_LWadj.LAT];
LON = [ITRDB_LWadj.LON];
plotm(LAT, LON, 'k^', 'MarkerSize',3, 'MarkerFaceColor',[208,209,230]/255);
axis image;

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','ITRDB_EW-LW.tif')
close all;

