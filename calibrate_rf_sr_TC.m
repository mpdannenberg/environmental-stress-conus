%% Run RF models for each ecoregion
% Relative stress index (Sr), predicted using topography and climate (TC)
% M Dannenberg, 15 Apr 2016

%% Parameters
n_min = 75;
syear = 1971;

% Number of vars for each environmental stress type
nTopo = 4;
nSoil = 0; 
nTemp = 8;
nWater = 8;

% Names of variables
TopoLabs = {'ELEV','SLOPE','UAA','TWI'};
SoilLabs = {''};
TempLabs = {'TMIN_{SON}','TMIN_{DJF}','TMIN_{MAM}','TMIN_{JJA}',...
	'TMAX_{SON}','TMAX_{DJF}','TMAX_{MAM}','TMAX_{JJA}'};
WaterLabs = {'WB_{SON}','WB_{DJF}','WB_{MAM}','WB_{JJA}',...
	'VPD_{SON}','VPD_{DJF}','VPD_{MAM}','VPD_{JJA}'};
AllLabs = [TopoLabs TempLabs WaterLabs];

p = length(AllLabs);

% month of previous year to start with
smonth = 9; % September


%% Load/process datasets
load TREESI_wEnvFactors_RSI;

EcoL1 = [TREESI.EcoL1];
EcoL2 = [TREESI.EcoL2];
EcoL3 = [TREESI.EcoL3];

L1list = sort(unique(EcoL1));
L1list = L1list(L1list >0);

L2list = sort(unique(EcoL2));
L2list = L2list(L2list >0);

L3list = sort(unique(EcoL3));
L3list = L3list(L3list >0);



%% Level 1 models
tic
TopoImpEcoL1 = NaN(length(L1list), nTopo);
ClimImpEcoL1 = NaN(length(L1list), nTemp+nWater);

L1_stats = NaN(length(L1list), 4); % ecoregion, rmse, r2, Spearman's rho

for i = L1list
    rowi = find(L1list == i);
    treesi_sub = TREESI(EcoL1==i);
    n = sum([treesi_sub.END] - (syear-1));
    
    nSites = sum(EcoL1==i);
    
    if n>=n_min
        Y = NaN(n, 1);
        X = NaN(n, (nTopo+nSoil+nTemp+nWater));
        
        % Add site-level data to tables
        idx = 1;
        for s = 1:nSites
            eyear = treesi_sub(s).END;
            yr = treesi_sub(s).START:treesi_sub(s).END;
            inds = find(yr>=syear);
            nyrs = length(inds);
            
            Y(idx:(idx+nyrs-1)) = treesi_sub(s).TREESI(inds);
            
            % Topography
            X(idx:(idx+nyrs-1), 1) = treesi_sub(s).ELEV;
            X(idx:(idx+nyrs-1), 2) = treesi_sub(s).SLOPE;
            X(idx:(idx+nyrs-1), 3) = treesi_sub(s).UAA;
            X(idx:(idx+nyrs-1), 4) = treesi_sub(s).TWI;
            
            % Minimum Temperature
            pind = find(treesi_sub(s).PRISM.YEAR >= syear-1 & treesi_sub(s).PRISM.YEAR <= eyear-1);
            cind = find(treesi_sub(s).PRISM.YEAR >= syear & treesi_sub(s).PRISM.YEAR <= eyear);
            temp = NaN(nyrs, 12);
            temp(:, 1:(12-smonth+1)) = treesi_sub(s).PRISM.TMIN(pind, smonth:12);
            temp(:, (12-smonth+2):12) = treesi_sub(s).PRISM.TMIN(cind, 1:(smonth-1));
            temp2 = NaN(nyrs, 4);
            for mnth = 3:3:12
                temp2(:, mnth/3) = mean(temp(:, (mnth-2):mnth), 2);
            end
            X(idx:(idx+nyrs-1), (nTopo+nSoil+1):(nTopo+nSoil+4)) = temp2;
            
            % Maximum Temperature
            temp = NaN(nyrs, 12);
            temp(:, 1:(12-smonth+1)) = treesi_sub(s).PRISM.TMAX(pind, smonth:12);
            temp(:, (12-smonth+2):12) = treesi_sub(s).PRISM.TMAX(cind, 1:(smonth-1));
            temp2 = NaN(nyrs, 4);
            for mnth = 3:3:12
                temp2(:, mnth/3) = mean(temp(:, (mnth-2):mnth), 2);
            end
            X(idx:(idx+nyrs-1), (nTopo+nSoil+5):(nTopo+nSoil+8)) = temp2;
            
            % Water balance (P-PET)
            water = treesi_sub(s).PRISM.PPT - treesi_sub(s).PRISM.PET;
            temp = NaN(nyrs, 12);
            temp(:, 1:(12-smonth+1)) = water(pind, smonth:12);
            temp(:, (12-smonth+2):12) = water(cind, 1:(smonth-1));
            temp2 = NaN(nyrs, 4);
            for mnth = 3:3:12
                temp2(:, mnth/3) = sum(temp(:, (mnth-2):mnth), 2);
            end
            X(idx:(idx+nyrs-1), (nTopo+nSoil+nTemp+1):(nTopo+nSoil+nTemp+4)) = temp2;
            
            % VPD
            temp = NaN(nyrs, 12);
            temp(:, 1:(12-smonth+1)) = treesi_sub(s).PRISM.VPD(pind, smonth:12);
            temp(:, (12-smonth+2):12) = treesi_sub(s).PRISM.VPD(cind, 1:(smonth-1));
            temp2 = NaN(nyrs, 4);
            for mnth = 3:3:12
                temp2(:, mnth/3) = mean(temp(:, (mnth-2):mnth), 2);
            end
            X(idx:(idx+nyrs-1), (nTopo+nSoil+nTemp+5):(nTopo+nSoil+nTemp+8)) = temp2;
            
            idx = idx+nyrs;
        end
        
        cd RF_v6_RSI_TC;
        % Train random forest
		B = TreeBagger(300, X, Y, 'Method','regression',...
            'OOBVarImp','on');
        Err = oobError(B);
        plot(sqrt(Err));
        set(gca, 'FontName','Times New Roman');
        xlabel('Number of grown trees', 'FontName','Times New Roman');
        ylabel('OOB Root Mean Squared Error', 'FontName','Times New Roman');
        title(['Ecoregion ', num2str(i)], 'FontName','Times New Roman');
        fname = ['OOBError_EcoL1_', num2str(i), '.tif'];
        set(gcf,'PaperPositionMode','auto')
        print('-dtiff','-f1','-r600',fname)
        close all;

        
        [imp, impidx] = sort(B.OOBPermutedVarDeltaError, 'descend');
        plot(1:length(imp), imp, 'ko-', 'LineWidth',1.5);
        set(gcf, 'Position',[335 319 1001 463]);
        set(gca,'XTick',1:p, 'XTickLabel',AllLabs(impidx), 'XLim',[0.5 p+0.5], 'FontName','Times New Roman');
        ax = gca;
        ax.XTickLabelRotation = 90;
        grid on;
        ylabel('Mean Decrease in Accuracy', 'FontName','Times New Roman');
        title(['Ecoregion ', num2str(i)], 'FontName','Times New Roman');
        fname = ['VarImportance_EcoL1_', num2str(i), '.tif'];
        set(gcf,'PaperPositionMode','auto')
        print('-dtiff','-f1','-r600',fname)
        close all;
        
        yhat = oobPredict(B);
        plot(B.Y, yhat, 'ko');
        set(gcf, 'Position', [520 293 560 505]);
        set(gca, 'XLim',[0 2], 'YLim',[0 2], 'FontName','Times New Roman');
        xlabel('Observed \itf(E)', 'FontName','Times New Roman');
        ylabel('Predicted \itf(E)', 'FontName','Times New Roman');
        title(['Ecoregion ', num2str(i)], 'FontName','Times New Roman')
        r2 = corr(B.Y, yhat, 'rows','complete')^2;
        text(0.1, 1.9, ['r^{2} = ', num2str(round(r2, 2))], 'FontName','Times New Roman');
        hold on;
        plot([0 2], [0 2], '-', 'Color',[0.6 0.6 0.6]);
        fname = ['ScatterPlot_EcoL1_', num2str(i), '.tif'];
        set(gcf,'PaperPositionMode','auto')
        print('-dtiff','-f1','-r600',fname)
        close all;

        cd ..
        
        TopoImpEcoL1(L1list==i, :) = B.OOBPermutedVarDeltaError(1:nTopo);
        ClimImpEcoL1(L1list==i, :) = B.OOBPermutedVarDeltaError((nTopo+nSoil+1):(nTopo+nSoil+nTemp+nWater));
        
        eval(['B_',num2str(i), '=B;']);
        
        L1_stats(rowi, 1) = i;
        L1_stats(rowi, 2) = sqrt(Err(end));
        L1_stats(rowi, 3) = r2;
        L1_stats(rowi, 4) = corr(B.Y, yhat, 'type','Spearman', 'rows','complete');
        
    end
    
    
    
end

toc


%% Level 2 models
tic
TopoImpEcoL2 = NaN(length(L2list), nTopo);
ClimImpEcoL2 = NaN(length(L2list), nTemp+nWater);

L2_stats = NaN(length(L2list), 4); % ecoregion, rmse, r2, Spearman's rho

for i = L2list
    rowi = find(L2list == i);
    treesi_sub = TREESI(EcoL2==i);
    n = sum([treesi_sub.END] - (syear-1));
    
    nSites = sum(EcoL2==i);
    
    if n>=n_min
        Y = NaN(n, 1);
        X = NaN(n, (nTopo+nSoil+nTemp+nWater));
        
        % Add site-level data to tables
        idx = 1;
        for s = 1:nSites
            eyear = treesi_sub(s).END;
            yr = treesi_sub(s).START:treesi_sub(s).END;
            inds = find(yr>=syear);
            nyrs = length(inds);
            
            Y(idx:(idx+nyrs-1)) = treesi_sub(s).TREESI(inds);
            
            % Topography
            X(idx:(idx+nyrs-1), 1) = treesi_sub(s).ELEV;
            X(idx:(idx+nyrs-1), 2) = treesi_sub(s).SLOPE;
            X(idx:(idx+nyrs-1), 3) = treesi_sub(s).UAA;
            X(idx:(idx+nyrs-1), 4) = treesi_sub(s).TWI;
            
            % Minimum Temperature
            pind = find(treesi_sub(s).PRISM.YEAR >= syear-1 & treesi_sub(s).PRISM.YEAR <= eyear-1);
            cind = find(treesi_sub(s).PRISM.YEAR >= syear & treesi_sub(s).PRISM.YEAR <= eyear);
            temp = NaN(nyrs, 12);
            temp(:, 1:(12-smonth+1)) = treesi_sub(s).PRISM.TMIN(pind, smonth:12);
            temp(:, (12-smonth+2):12) = treesi_sub(s).PRISM.TMIN(cind, 1:(smonth-1));
            temp2 = NaN(nyrs, 4);
            for mnth = 3:3:12
                temp2(:, mnth/3) = mean(temp(:, (mnth-2):mnth), 2);
            end
            X(idx:(idx+nyrs-1), (nTopo+nSoil+1):(nTopo+nSoil+4)) = temp2;
            
            % Maximum Temperature
            temp = NaN(nyrs, 12);
            temp(:, 1:(12-smonth+1)) = treesi_sub(s).PRISM.TMAX(pind, smonth:12);
            temp(:, (12-smonth+2):12) = treesi_sub(s).PRISM.TMAX(cind, 1:(smonth-1));
            temp2 = NaN(nyrs, 4);
            for mnth = 3:3:12
                temp2(:, mnth/3) = mean(temp(:, (mnth-2):mnth), 2);
            end
            X(idx:(idx+nyrs-1), (nTopo+nSoil+5):(nTopo+nSoil+8)) = temp2;
            
            % Water balance (P-PET)
            water = treesi_sub(s).PRISM.PPT - treesi_sub(s).PRISM.PET;
            temp = NaN(nyrs, 12);
            temp(:, 1:(12-smonth+1)) = water(pind, smonth:12);
            temp(:, (12-smonth+2):12) = water(cind, 1:(smonth-1));
            temp2 = NaN(nyrs, 4);
            for mnth = 3:3:12
                temp2(:, mnth/3) = sum(temp(:, (mnth-2):mnth), 2);
            end
            X(idx:(idx+nyrs-1), (nTopo+nSoil+nTemp+1):(nTopo+nSoil+nTemp+4)) = temp2;
            
            % VPD
            temp = NaN(nyrs, 12);
            temp(:, 1:(12-smonth+1)) = treesi_sub(s).PRISM.VPD(pind, smonth:12);
            temp(:, (12-smonth+2):12) = treesi_sub(s).PRISM.VPD(cind, 1:(smonth-1));
            temp2 = NaN(nyrs, 4);
            for mnth = 3:3:12
                temp2(:, mnth/3) = mean(temp(:, (mnth-2):mnth), 2);
            end
            X(idx:(idx+nyrs-1), (nTopo+nSoil+nTemp+5):(nTopo+nSoil+nTemp+8)) = temp2;
            
            idx = idx+nyrs;
        end
        
        cd RF_v6_RSI_TC;
        % Train random forest
		B = TreeBagger(300, X, Y, 'Method','regression',...
            'OOBVarImp','on');
        Err = oobError(B);
        plot(sqrt(Err));
        set(gca, 'FontName','Times New Roman');
        xlabel('Number of grown trees', 'FontName','Times New Roman');
        ylabel('OOB Root Mean Squared Error', 'FontName','Times New Roman');
        title(['Ecoregion ', num2str(i)], 'FontName','Times New Roman');
        fname = ['OOBError_EcoL2_', num2str(i), '.tif'];
        set(gcf,'PaperPositionMode','auto')
        print('-dtiff','-f1','-r600',fname)
        close all;

        
        [imp, impidx] = sort(B.OOBPermutedVarDeltaError, 'descend');
        plot(1:length(imp), imp, 'ko-', 'LineWidth',1.5);
        set(gcf, 'Position',[335 319 1001 463]);
        set(gca,'XTick',1:p, 'XTickLabel',AllLabs(impidx), 'XLim',[0.5 p+0.5], 'FontName','Times New Roman');
        ax = gca;
        ax.XTickLabelRotation = 90;
        grid on;
        ylabel('Mean Decrease in Accuracy', 'FontName','Times New Roman');
        title(['Ecoregion ', num2str(i)], 'FontName','Times New Roman');
        fname = ['VarImportance_EcoL2_', num2str(i), '.tif'];
        set(gcf,'PaperPositionMode','auto')
        print('-dtiff','-f1','-r600',fname)
        close all;
        
        yhat = oobPredict(B);
        plot(B.Y, yhat, 'ko');
        set(gcf, 'Position', [520 293 560 505]);
        set(gca, 'XLim',[0 2], 'YLim',[0 2], 'FontName','Times New Roman');
        xlabel('Observed \itf(E)', 'FontName','Times New Roman');
        ylabel('Predicted \itf(E)', 'FontName','Times New Roman');
        title(['Ecoregion ', num2str(i)], 'FontName','Times New Roman')
        r2 = corr(B.Y, yhat, 'rows','complete')^2;
        text(0.1, 1.9, ['r^{2} = ', num2str(round(r2, 2))], 'FontName','Times New Roman');
        hold on;
        plot([0 2], [0 2], '-', 'Color',[0.6 0.6 0.6]);
        fname = ['ScatterPlot_EcoL2_', num2str(i), '.tif'];
        set(gcf,'PaperPositionMode','auto')
        print('-dtiff','-f1','-r600',fname)
        close all;

        cd ..
        
        TopoImpEcoL2(L2list==i, :) = B.OOBPermutedVarDeltaError(1:nTopo);
        ClimImpEcoL2(L2list==i, :) = B.OOBPermutedVarDeltaError((nTopo+nSoil+1):(nTopo+nSoil+nTemp+nWater));
        
        eval(['B_',num2str(i), '=B;']);
        
        L2_stats(rowi, 1) = i;
        L2_stats(rowi, 2) = sqrt(Err(end));
        L2_stats(rowi, 3) = r2;
        L2_stats(rowi, 4) = corr(B.Y, yhat, 'type','Spearman', 'rows','complete');
        
    end
    
    
    
end

toc


%% Level 3 models
tic
TopoImpEcoL3 = NaN(length(L3list), nTopo);
ClimImpEcoL3 = NaN(length(L3list), nTemp+nWater);

L3_stats = NaN(length(L3list), 4); % ecoregion, rmse, r2, Spearman's rho

for i = L3list
    rowi = find(L3list == i);
    treesi_sub = TREESI(EcoL3==i);
    n = sum([treesi_sub.END] - (syear-1));
    
    nSites = sum(EcoL3==i);
    
    if n>=n_min
        Y = NaN(n, 1);
        X = NaN(n, (nTopo+nSoil+nTemp+nWater));
        
        % Add site-level data to tables
        idx = 1;
        for s = 1:nSites
            eyear = treesi_sub(s).END;
            yr = treesi_sub(s).START:treesi_sub(s).END;
            inds = find(yr>=syear);
            nyrs = length(inds);
            
            Y(idx:(idx+nyrs-1)) = treesi_sub(s).TREESI(inds);
            
            % Topography
            X(idx:(idx+nyrs-1), 1) = treesi_sub(s).ELEV;
            X(idx:(idx+nyrs-1), 2) = treesi_sub(s).SLOPE;
            X(idx:(idx+nyrs-1), 3) = treesi_sub(s).UAA;
            X(idx:(idx+nyrs-1), 4) = treesi_sub(s).TWI;
            
            % Minimum Temperature
            pind = find(treesi_sub(s).PRISM.YEAR >= syear-1 & treesi_sub(s).PRISM.YEAR <= eyear-1);
            cind = find(treesi_sub(s).PRISM.YEAR >= syear & treesi_sub(s).PRISM.YEAR <= eyear);
            temp = NaN(nyrs, 12);
            temp(:, 1:(12-smonth+1)) = treesi_sub(s).PRISM.TMIN(pind, smonth:12);
            temp(:, (12-smonth+2):12) = treesi_sub(s).PRISM.TMIN(cind, 1:(smonth-1));
            temp2 = NaN(nyrs, 4);
            for mnth = 3:3:12
                temp2(:, mnth/3) = mean(temp(:, (mnth-2):mnth), 2);
            end
            X(idx:(idx+nyrs-1), (nTopo+nSoil+1):(nTopo+nSoil+4)) = temp2;
            
            % Maximum Temperature
            temp = NaN(nyrs, 12);
            temp(:, 1:(12-smonth+1)) = treesi_sub(s).PRISM.TMAX(pind, smonth:12);
            temp(:, (12-smonth+2):12) = treesi_sub(s).PRISM.TMAX(cind, 1:(smonth-1));
            temp2 = NaN(nyrs, 4);
            for mnth = 3:3:12
                temp2(:, mnth/3) = mean(temp(:, (mnth-2):mnth), 2);
            end
            X(idx:(idx+nyrs-1), (nTopo+nSoil+5):(nTopo+nSoil+8)) = temp2;
            
            % Water balance (P-PET)
            water = treesi_sub(s).PRISM.PPT - treesi_sub(s).PRISM.PET;
            temp = NaN(nyrs, 12);
            temp(:, 1:(12-smonth+1)) = water(pind, smonth:12);
            temp(:, (12-smonth+2):12) = water(cind, 1:(smonth-1));
            temp2 = NaN(nyrs, 4);
            for mnth = 3:3:12
                temp2(:, mnth/3) = sum(temp(:, (mnth-2):mnth), 2);
            end
            X(idx:(idx+nyrs-1), (nTopo+nSoil+nTemp+1):(nTopo+nSoil+nTemp+4)) = temp2;
            
            % VPD
            temp = NaN(nyrs, 12);
            temp(:, 1:(12-smonth+1)) = treesi_sub(s).PRISM.VPD(pind, smonth:12);
            temp(:, (12-smonth+2):12) = treesi_sub(s).PRISM.VPD(cind, 1:(smonth-1));
            temp2 = NaN(nyrs, 4);
            for mnth = 3:3:12
                temp2(:, mnth/3) = mean(temp(:, (mnth-2):mnth), 2);
            end
            X(idx:(idx+nyrs-1), (nTopo+nSoil+nTemp+5):(nTopo+nSoil+nTemp+8)) = temp2;
            
            idx = idx+nyrs;
        end
        
        cd RF_v6_RSI_TC;
        % Train random forest
		B = TreeBagger(300, X, Y, 'Method','regression',...
            'OOBVarImp','on');
        Err = oobError(B);
        plot(sqrt(Err));
        set(gca, 'FontName','Times New Roman');
        xlabel('Number of grown trees', 'FontName','Times New Roman');
        ylabel('OOB Root Mean Squared Error', 'FontName','Times New Roman');
        title(['Ecoregion ', num2str(i)], 'FontName','Times New Roman');
        fname = ['OOBError_EcoL3_', num2str(i), '.tif'];
        set(gcf,'PaperPositionMode','auto')
        print('-dtiff','-f1','-r600',fname)
        close all;

        
        [imp, impidx] = sort(B.OOBPermutedVarDeltaError, 'descend');
        plot(1:length(imp), imp, 'ko-', 'LineWidth',1.5);
        set(gcf, 'Position',[335 319 1001 463]);
        set(gca,'XTick',1:p, 'XTickLabel',AllLabs(impidx), 'XLim',[0.5 p+0.5], 'FontName','Times New Roman');
        ax = gca;
        ax.XTickLabelRotation = 90;
        grid on;
        ylabel('Mean Decrease in Accuracy', 'FontName','Times New Roman');
        title(['Ecoregion ', num2str(i)], 'FontName','Times New Roman');
        fname = ['VarImportance_EcoL3_', num2str(i), '.tif'];
        set(gcf,'PaperPositionMode','auto')
        print('-dtiff','-f1','-r600',fname)
        close all;
        
        yhat = oobPredict(B);
        plot(B.Y, yhat, 'ko');
        set(gcf, 'Position', [520 293 560 505]);
        set(gca, 'XLim',[0 2], 'YLim',[0 2], 'FontName','Times New Roman');
        xlabel('Observed \itf(E)', 'FontName','Times New Roman');
        ylabel('Predicted \itf(E)', 'FontName','Times New Roman');
        title(['Ecoregion ', num2str(i)], 'FontName','Times New Roman')
        r2 = corr(B.Y, yhat, 'rows','complete')^2;
        text(0.1, 1.9, ['r^{2} = ', num2str(round(r2, 2))], 'FontName','Times New Roman');
        hold on;
        plot([0 2], [0 2], '-', 'Color',[0.6 0.6 0.6]);
        fname = ['ScatterPlot_EcoL3_', num2str(i), '.tif'];
        set(gcf,'PaperPositionMode','auto')
        print('-dtiff','-f1','-r600',fname)
        close all;

        cd ..
        
        TopoImpEcoL3(L3list==i, :) = B.OOBPermutedVarDeltaError(1:nTopo);
        ClimImpEcoL3(L3list==i, :) = B.OOBPermutedVarDeltaError((nTopo+nSoil+1):(nTopo+nSoil+nTemp+nWater));
        
        eval(['B_',num2str(i), '=B;']);
        
        L3_stats(rowi, 1) = i;
        L3_stats(rowi, 2) = sqrt(Err(end));
        L3_stats(rowi, 3) = r2;
        L3_stats(rowi, 4) = corr(B.Y, yhat, 'type','Spearman', 'rows','complete');
        
    end
    
    
    
end

toc

%% Importance boxplots
cd RF_v6_RSI_TC;

% Topography importance
h = figure('Color','w');
set(h, 'Position',[520 527 846 271]);

ha = tight_subplot(1,3,0.02, [0.2 0.1],[0.1 0.05]);
axes(ha(1));
boxplot(TopoImpEcoL1, 'PlotStyle','compact', 'Colors',[0.4 0.4 0.4],...
    'Orientation','Horizontal', 'Jitter',0);
ax = gca;
set(ax, 'FontName','Times New Roman', 'YTick', 1:nTopo,'XLim',[-0.3 1.5],... 
    'YTickLabel',TopoLabs);
n = sum(isfinite(TopoImpEcoL1(:,1)));
title(['Level 1 (','{\itn} = ',num2str(n), ')'], 'FontName','Times New Roman');
grid on;

axes(ha(2));
boxplot(TopoImpEcoL2, 'PlotStyle','compact', 'Colors',[0.4 0.4 0.4],...
    'Orientation','Horizontal', 'Jitter',0);
ax = gca;
set(ax, 'FontName','Times New Roman', 'YTick', 1:nTopo,'XLim',[-0.3 1.5],...
    'YTickLabel',{'','','',''});
n = sum(isfinite(TopoImpEcoL2(:,1)));
title(['Level 2 (','{\itn} = ',num2str(n), ')'], 'FontName','Times New Roman');
text(0.5,-0.3,'Mean Decrease in Accuracy when Variable is Excluded',...
    'FontName','Times New Roman', 'HorizontalAlignment','center', 'FontSize',14);
grid on;

axes(ha(3));
boxplot(TopoImpEcoL3, 'PlotStyle','compact', 'Colors',[0.4 0.4 0.4],...
    'Orientation','Horizontal', 'Jitter',0);
ax = gca;
set(ax, 'FontName','Times New Roman', 'YTick', 1:nTopo,'XLim',[-0.3 1.5],...
    'YTickLabel',{'','','',''});
n = sum(isfinite(TopoImpEcoL3(:,1)));
title(['Level 3 (','{\itn} = ',num2str(n), ')'], 'FontName','Times New Roman');
grid on;

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r600','TopoImportanceBox.tif')

% Climate importance
h = figure('Color','w');
set(h, 'Position',[369 287 997 511]);

ha = tight_subplot(1,3,0.02, [0.1 0.05],[0.1 0.05]);
axes(ha(1));
boxplot(ClimImpEcoL1, 'PlotStyle','compact', 'Colors',[0.4 0.4 0.4],...
    'Orientation','Horizontal', 'Jitter',0);
ax = gca;
set(ax, 'FontName','Times New Roman', 'YTick', 1:(nTemp+nWater),'XLim',[-0.4 1.45],... 
    'YTickLabel',[TempLabs WaterLabs]);
n = sum(isfinite(ClimImpEcoL1(:,1)));
title(['Level 1 (','{\itn} = ',num2str(n), ')'], 'FontName','Times New Roman');
grid on;

axes(ha(2));
boxplot(ClimImpEcoL2, 'PlotStyle','compact', 'Colors',[0.4 0.4 0.4],...
    'Orientation','Horizontal', 'Jitter',0);
ax = gca;
set(ax, 'FontName','Times New Roman', 'YTick', 1:(nTemp+nWater),'XLim',[-0.4 1.45],...
    'YTickLabel',{'','','',''});
n = sum(isfinite(ClimImpEcoL2(:,1)));
title(['Level 2 (','{\itn} = ',num2str(n), ')'], 'FontName','Times New Roman');
text(0.5,-0.9,'Mean Decrease in Accuracy when Variable is Excluded',...
    'FontName','Times New Roman', 'HorizontalAlignment','center', 'FontSize',14);
grid on;

axes(ha(3));
boxplot(ClimImpEcoL3, 'PlotStyle','compact', 'Colors',[0.4 0.4 0.4],...
    'Orientation','Horizontal', 'Jitter',0);
ax = gca;
set(ax, 'FontName','Times New Roman', 'YTick', 1:(nTemp+nWater),'XLim',[-0.4 1.45],...
    'YTickLabel',{'','','',''});
n = sum(isfinite(ClimImpEcoL3(:,1)));
title(['Level 3 (','{\itn} = ',num2str(n), ')'], 'FontName','Times New Roman');
grid on;

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f2','-r600','ClimImportanceBox.tif')

close all;


%% Validation Maps
cd ..;
load conus_mask;

cd EcoRegions;
[EcoL1, R] = geotiffread('us_EcoL1_4km.tif');
[EcoL2] = geotiffread('us_EcoL2_4km.tif');
[EcoL3] = geotiffread('us_EcoL3_4km.tif');

RF_RMSE = NaN(size(EcoL1));
RF_R2 = NaN(size(EcoL1));
RF_RHO = NaN(size(EcoL1));

% Level 1
for i = 1:length(L1list)
    
    RF_RMSE(EcoL1 == L1_stats(i, 1)) = L1_stats(i, 2);
    RF_R2(EcoL1 == L1_stats(i, 1)) = L1_stats(i, 3);
    RF_RHO(EcoL1 == L1_stats(i, 1)) = L1_stats(i, 4);
    
end

% Level 2
for i = 1:length(L2list)
    
    RF_RMSE(EcoL2 == L2_stats(i, 1)) = L2_stats(i, 2);
    RF_R2(EcoL2 == L2_stats(i, 1)) = L2_stats(i, 3);
    RF_RHO(EcoL2 == L2_stats(i, 1)) = L2_stats(i, 4);
    
end

% Level 3
for i = 1:length(L3list)
    
    RF_RMSE(EcoL3 == L3_stats(i, 1)) = L3_stats(i, 2);
    RF_R2(EcoL3 == L3_stats(i, 1)) = L3_stats(i, 3);
    RF_RHO(EcoL3 == L3_stats(i, 1)) = L3_stats(i, 4);
    
end

cd ../PRISM

load prism_latlon;
cd ..;

states = shaperead('usastatehi','UseGeoCoords', true);
latlim = [25 50];
lonlim = [-126 -65];

h=figure('Color','w');
set(h, 'Position',[9 49 556 768]);

% RMSE
subplot(3, 1, 1)
ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',10,'MLineLocation',20,'MeridianLabel','on',...
        'ParallelLabel','on','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Times New Roman',...
        'FLineWidth',1);
h=surfm(PRISMlat, PRISMlon, RF_RMSE.*CONUS);
geoshow(states,'FaceColor','none','LineWidth',0.4)
caxis([0 0.25])
colormap(flipud(gray(10)))
axis off;
axis image;
freezeColors
h = colorbar('eastoutside');
pos = get(h, 'Position');
pos(1) = pos(1)+0.02;
set(h, 'Position',pos, 'FontName','Times New Roman', 'YTick',0:0.05:0.25);
set(ax, 'Position', [0.05 0.7093 0.7750 0.2157])
ylabel(h, 'RMSE', 'FontName','Times New Roman', 'FontSize',14);
text(-0.58,0.9, 'a)', 'FontName','Times New Roman', 'FontSize',18);

% R^2
subplot(3, 1, 2)
ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',10,'MLineLocation',20,'MeridianLabel','on',...
        'ParallelLabel','on','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Times New Roman',...
        'FLineWidth',1);
h=surfm(PRISMlat, PRISMlon, RF_R2.*CONUS);
geoshow(states,'FaceColor','none','LineWidth',0.4)
caxis([0 1])
colormap(flipud(gray(10)))
axis off;
axis image;
freezeColors
h = colorbar('eastoutside');
pos = get(h, 'Position');
pos(1) = pos(1)+0.02;
set(h, 'Position',pos, 'FontName','Times New Roman', 'YTick',0:0.2:1);
set(ax, 'Position', [0.05 0.4096 0.7750 0.2157])
ylabel(h, 'r^{2}', 'FontName','Times New Roman', 'FontSize',14);
text(-0.58,0.9, 'b)', 'FontName','Times New Roman', 'FontSize',18);

% Spearman's rho
subplot(3, 1, 3)
ax = axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',10,'MLineLocation',20,'MeridianLabel','on',...
        'ParallelLabel','on','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Times New Roman',...
        'FLineWidth',1);
h=surfm(PRISMlat, PRISMlon, RF_RHO.*CONUS);
geoshow(states,'FaceColor','none','LineWidth',0.4)
caxis([0 1])
colormap(flipud(gray(10)))
axis off;
axis image;
freezeColors
h = colorbar('eastoutside');
pos = get(h, 'Position');
pos(1) = pos(1)+0.02;
set(h, 'Position',pos, 'FontName','Times New Roman', 'YTick',0:0.2:1);
set(ax, 'Position', [0.05 0.1100 0.7750 0.2157])
ylabel(h, 'Spearman''s \rho', 'FontName','Times New Roman', 'FontSize',14);
text(-0.58,0.9, 'c)', 'FontName','Times New Roman', 'FontSize',18);

cd RF_v6_RSI_TC;
set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r600','RF_v6_RSI_TC_ValStats.tif')

save('RF_v6_RSI_TC_stats.mat','L1_stats','L2_stats','L3_stats', '-v7.3');

cd ..
close all;

