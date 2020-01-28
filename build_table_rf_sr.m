% Build table of RF validation statistics for relative stress index

% Full Model (TSC)
cd RF_v6_RSI_TSC;
load RF_v6_RSI_TSC_stats.mat;

L1_stats = L1_stats(~any(isnan(L1_stats),2),:);
L2_stats = L2_stats(~any(isnan(L2_stats),2),:);
L3_stats = L3_stats(~any(isnan(L3_stats),2),:);

A = [L1_stats;L2_stats;L3_stats];

T = table;
T.Ecoregion = A(:, 1);
T.Sites = A(:, 6);
T.n = A(:,5);
T.RMSE = A(:,2);
T.r2 = A(:,3);
T.rho = A(:,4);

clearvars -except T;

% TC
cd ../RF_v6_RSI_TC;
load RF_v6_RSI_TC_stats.mat;

L1_stats = L1_stats(~any(isnan(L1_stats),2),:);
L2_stats = L2_stats(~any(isnan(L2_stats),2),:);
L3_stats = L3_stats(~any(isnan(L3_stats),2),:);

A = [L1_stats;L2_stats;L3_stats];

T.dr2_TC = A(:, 3) - T.r2;

clearvars -except T;

% SC
cd ../RF_v6_RSI_SC;
load RF_v6_RSI_SC_stats.mat;

L1_stats = L1_stats(~any(isnan(L1_stats),2),:);
L2_stats = L2_stats(~any(isnan(L2_stats),2),:);
L3_stats = L3_stats(~any(isnan(L3_stats),2),:);

A = [L1_stats;L2_stats;L3_stats];

T.dr2_SC = A(:, 3) - T.r2;


clearvars -except T;

% C
cd ../RF_v6_RSI_C;
load RF_v6_RSI_C_stats.mat;

L1_stats = L1_stats(~any(isnan(L1_stats),2),:);
L2_stats = L2_stats(~any(isnan(L2_stats),2),:);
L3_stats = L3_stats(~any(isnan(L3_stats),2),:);

A = [L1_stats;L2_stats;L3_stats];

T.dr2_C = A(:, 3) - T.r2;

cd ..


