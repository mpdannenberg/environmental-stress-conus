% Read stress index files produced from ITRDB rwls
% M Dannenberg 8/9/14
% modified 9/12/15

tic

eyear = 2013; % Don't use anything after this (no climate data)

[ninfo, tinfo, ~] = xlsread('ITRDB_NAm_091115.xlsx','B:K');

cd StressIndex; % Folder with stress index chronologies from calculate_stress_index_itrdb.R
fnames = glob('*.csv');

SITE = tinfo(2:end, 1);
SPECIES = tinfo(2:end, 4);
START = ninfo(:, 1);
END = ninfo(:, 2);
LAT = ninfo(:, 4);
LON = ninfo(:, 5);
ELEV = ninfo(:, 6);
clear ninfo tinfo;

for i = 1:length(fnames)
    sname = regexp(fnames{i}, '\_', 'split');
    TF = strcmpi(sname{1}, SITE);
    idx = find(TF);
    
    % Basic Info
    TREESI(i).SITE = SITE{idx};
    TREESI(i).SPECIES = SPECIES{idx};
    TREESI(i).START = START(idx);
    TREESI(i).END = END(idx);
    TREESI(i).ELEV = ELEV(idx);
    
    try
        % Degree Minute --> Decimal degrees
        dsign = sign(LAT(idx));
        degree = floor(abs(LAT(idx)));
        minute = round(100*(abs(LAT(idx))-degree));
        TREESI(i).LAT = dsign * dms2degrees([degree minute 0]);

        dsign = sign(LON(idx));
        degree = floor(abs(LON(idx))); 
        minute = round(100*(abs(LON(idx))-degree));
        TREESI(i).LON = dsign * dms2degrees([degree minute 0]);
        clear dsign degree minute;
    catch
        TREESI(i).LAT = LAT(idx);
        TREESI(i).LON = LON(idx);
    end
    
    % Read CSV with RWI
    out = read_mixed_csv(fnames{i}, ',');
    trsi = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,2)));
    samp = cell2mat(cellfun(@(s) {str2double(s)}, out(2:end,3)));
    year = str2double(strrep(out(2:end, 1), '"', ''));
    
    TREESI(i).TREESI = trsi(year<=eyear);
    TREESI(i).SAMPLE_DEPTH = samp(year<=eyear);
    TREESI(i).YEAR = year(year<=eyear);
    TREESI(i).START = min(year);
    TREESI(i).END = min(max(year), 2013);
    
    clear out trsi samp year;
    
    clear sname TF idx;
    
end

clear ELEV END LAT LOCATION LON NAME PI SITE SPECIES START fnames i;

cd ..

save('./data/TREESI.mat', 'TREESI');

toc


