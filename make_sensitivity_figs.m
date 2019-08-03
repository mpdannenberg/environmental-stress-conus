% New supplemental figures for Global Ecol & Biogeogr paper
clr = wesanderson('fantasticfox1');

%% Bias in DBH measurements (Fig. A12 of diss.)

pipo = csvread('./data/DBH_Comparison_PIPO.csv', 1, 2);
psme = csvread('./data/DBH_Comparison_PSME.csv', 1, 2);
acru = csvread('./data/DBH_Comparison_ACRU.csv', 1, 2);
quru = csvread('./data/DBH_Comparison_QURU.csv', 1, 2);
tsca = csvread('./data/DBH_Comparison_TSCA.csv', 1, 2);

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 5 6];

subplot(3,1,[1 2])
plot(pipo(:, 2), pipo(:, 1), 'o', 'Color',clr(1, :), 'LineWidth',...
    1.2, 'MarkerSize',5);
set(gca, 'XLim',[0 150], 'YLim',[0 150], 'TickDir','out', 'TickLength',[0.025 0.05], 'FontSize',9);
box off;
hold on;
plot([0 150], [0 150], '-', 'Color',[0.6 0.6 0.6]);
plot(psme(:, 2), psme(:, 1), 'o', 'Color',clr(2, :), 'LineWidth',...
    1.2, 'MarkerSize',5);
plot(quru(:, 2), quru(:, 1), 'o', 'Color',clr(4, :), 'LineWidth',...
    1.2, 'MarkerSize',5);
plot(acru(:, 2), acru(:, 1), 'o', 'Color',clr(3, :), 'LineWidth',...
    1.2, 'MarkerSize',5);
plot(tsca(:, 2), tsca(:, 1), 'o', 'Color',clr(5, :), 'LineWidth',...
    1.2, 'MarkerSize',5);
text(5, 142, '\itP. ponderosa', 'FontSize',11, 'Color',clr(1, :))
text(5, 132, '\itP. menziesii', 'FontSize',11, 'Color',clr(2, :))
text(5, 122, '\itA. rubrum', 'FontSize',11, 'Color',clr(3, :))
text(5, 112, '\itQ. rubra', 'FontSize',11, 'Color',clr(4, :))
text(5, 102, '\itT. canadensis', 'FontSize',11, 'Color',clr(5, :))
xlabel('Measured DBH (cm)', 'FontSize',12);
ylabel('Estimated DBH (cm)', 'FontSize',12);
set(gca, 'Position',[0.1300    0.43    0.7750    0.5154]);
text(-23, 150, 'A', 'FontSize',14, 'FontWeight','bold');

subplot(3, 1, 3)
Y1 = histcounts(pipo(:, 1)-pipo(:, 2), -100:5:20);
Y2 = histcounts(psme(:, 1)-psme(:, 2), -100:5:20);
Y3 = histcounts(acru(:, 1)-acru(:, 2), -100:5:20);
Y4 = histcounts(quru(:, 1)-quru(:, 2), -100:5:20);
Y5 = histcounts(tsca(:, 1)-tsca(:, 2), -100:5:20);
Y = [Y1' Y2' Y3' Y4' Y5'];
b = bar((-100+2.5):5:(20-2.5), Y, 1, 'stacked', 'EdgeColor',[0.4 0.4 0.4]);
b(1).FaceColor = clr(1, :);
b(2).FaceColor = clr(2, :);
b(3).FaceColor = clr(3, :);
b(4).FaceColor = clr(4, :);
b(5).FaceColor = clr(5, :);
set(gca, 'YLim',[0 320], 'TickDir','out', 'TickLength',[0.025 0.05], 'FontSize',9);
set(gca, 'Position',[0.1300    0.09    0.7750    0.2157]);
box off;
xlabel('DBH bias (cm)', 'FontSize',12);
ylabel('Count', 'FontSize',12);
text(-124, 320, 'B', 'FontSize',14, 'FontWeight','bold');

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/dbh-bias.tif')
close all;


%% Print stats:
% PIPO
fprintf('\n');
fprintf('PIPO:\n');
fprintf(['Mean bias: ', num2str(mean(pipo(:, 1)-pipo(:, 2))),'\n']);
fprintf(['Median bias: ', num2str(median(pipo(:, 1)-pipo(:, 2))),'\n']);
fprintf('\n');

% PSME
fprintf('PSME:\n');
fprintf(['Mean bias: ', num2str(mean(psme(:, 1)-psme(:, 2))),'\n']);
fprintf(['Median bias: ', num2str(median(psme(:, 1)-psme(:, 2))),'\n']);
fprintf('\n');

% ACRU
fprintf('\n');
fprintf('ACRU:\n');
fprintf(['Mean bias: ', num2str(mean(acru(:, 1)-acru(:, 2))),'\n']);
fprintf(['Median bias: ', num2str(median(acru(:, 1)-acru(:, 2))),'\n']);
fprintf('\n');

% QURU
fprintf('QURU:\n');
fprintf(['Mean bias: ', num2str(mean(quru(:, 1)-quru(:, 2))),'\n']);
fprintf(['Median bias: ', num2str(median(quru(:, 1)-quru(:, 2))),'\n']);
fprintf('\n');

% TSCA
fprintf('TSCA:\n');
fprintf(['Mean bias: ', num2str(mean(tsca(:, 1)-tsca(:, 2))),'\n']);
fprintf(['Median bias: ', num2str(median(tsca(:, 1)-tsca(:, 2))),'\n']);
fprintf('\n');


%% Sensitivity figures
nbins = 100; % Number of bins within which to calculate median and quantiles
quants = [0.025 0.975];

% 1-D*, 2-D, 3-Dopt*, 4-Dopt, 5-dD, 6-Dopt bias, 7-S*, 8-S, 9-S bias
pipo = csvread('./data/Sensitivity_PIPO.csv', 1, 1);
psme = csvread('./data/Sensitivity_PSME.csv', 1, 1);
acru = csvread('./data/Sensitivity_ACRU.csv', 1, 1);
quru = csvread('./data/Sensitivity_QURU.csv', 1, 1);
tsca = csvread('./data/Sensitivity_TSCA.csv', 1, 1);

% Plots: PIPO
Dstar = pipo(:,1);
Dopt_bias = pipo(:, 6);
S_bias = pipo(:, 9);
S_bias(~isfinite(S_bias)) = NaN;

PIPO.Dopt.y = NaN(1, nbins);
PIPO.Dopt.uci = NaN(1, nbins);
PIPO.Dopt.lci = NaN(1, nbins);
PIPO.S.y = NaN(1, nbins);
PIPO.S.uci = NaN(1, nbins);
PIPO.S.lci = NaN(1, nbins);

x = quantile(Dstar, 0:(1/nbins):1);
bins = (0+(0.5/nbins)):(1/nbins):(1-(0.5/nbins));
bins = quantile(Dstar, bins);
PIPO.Dopt.bins = bins;
PIPO.S.bins = bins;

for i = 1:nbins
    
    if i==1
        idx = Dstar>=x(i) & Dstar<=x(i+1);
    else
        idx = Dstar>x(i) & Dstar<=x(i+1);
    end
    Dopt_sub = Dopt_bias(idx);
    S_sub = S_bias(idx);
    
    PIPO.Dopt.y(i) = median(Dopt_sub);
    PIPO.S.y(i) = median(S_sub);
    
    cis = quantile(Dopt_sub, quants);
    PIPO.Dopt.lci(i) = cis(1);
    PIPO.Dopt.uci(i) = cis(2);
    
    cis = quantile(S_sub, quants);
    PIPO.S.lci(i) = cis(1);
    PIPO.S.uci(i) = cis(2);
    
end

% Plots: PSME
Dstar = psme(:,1);
Dopt_bias = psme(:, 6);
S_bias = psme(:, 9);
S_bias(~isfinite(S_bias)) = NaN;

PSME.Dopt.y = NaN(1, nbins);
PSME.Dopt.uci = NaN(1, nbins);
PSME.Dopt.lci = NaN(1, nbins);
PSME.S.y = NaN(1, nbins);
PSME.S.uci = NaN(1, nbins);
PSME.S.lci = NaN(1, nbins);

x = quantile(Dstar, 0:(1/nbins):1);
bins = (0+(0.5/nbins)):(1/nbins):(1-(0.5/nbins));
bins = quantile(Dstar, bins);
PSME.Dopt.bins = bins;
PSME.S.bins = bins;

for i = 1:nbins
    
    if i==1
        idx = Dstar>=x(i) & Dstar<=x(i+1);
    else
        idx = Dstar>x(i) & Dstar<=x(i+1);
    end
    Dopt_sub = Dopt_bias(idx);
    S_sub = S_bias(idx);
    
    PSME.Dopt.y(i) = median(Dopt_sub);
    PSME.S.y(i) = median(S_sub);
    
    cis = quantile(Dopt_sub, quants);
    PSME.Dopt.lci(i) = cis(1);
    PSME.Dopt.uci(i) = cis(2);
    
    cis = quantile(S_sub, quants);
    PSME.S.lci(i) = cis(1);
    PSME.S.uci(i) = cis(2);
    
end

% Plots: ACRU
Dstar = acru(:,1);
Dopt_bias = acru(:, 6);
S_bias = acru(:, 9);
S_bias(~isfinite(S_bias)) = NaN;

ACRU.Dopt.y = NaN(1, nbins);
ACRU.Dopt.uci = NaN(1, nbins);
ACRU.Dopt.lci = NaN(1, nbins);
ACRU.S.y = NaN(1, nbins);
ACRU.S.uci = NaN(1, nbins);
ACRU.S.lci = NaN(1, nbins);

x = quantile(Dstar, 0:(1/nbins):1);
bins = (0+(0.5/nbins)):(1/nbins):(1-(0.5/nbins));
bins = quantile(Dstar, bins);
ACRU.Dopt.bins = bins;
ACRU.S.bins = bins;

for i = 1:nbins
    
    if i==1
        idx = Dstar>=x(i) & Dstar<=x(i+1);
    else
        idx = Dstar>x(i) & Dstar<=x(i+1);
    end
    Dopt_sub = Dopt_bias(idx);
    S_sub = S_bias(idx);
    
    ACRU.Dopt.y(i) = median(Dopt_sub);
    ACRU.S.y(i) = median(S_sub);
    
    cis = quantile(Dopt_sub, quants);
    ACRU.Dopt.lci(i) = cis(1);
    ACRU.Dopt.uci(i) = cis(2);
    
    cis = quantile(S_sub, quants);
    ACRU.S.lci(i) = cis(1);
    ACRU.S.uci(i) = cis(2);
    
end

% Plots: QURU
Dstar = quru(:,1);
Dopt_bias = quru(:, 6);
S_bias = quru(:, 9);
S_bias(~isfinite(S_bias)) = NaN;

QURU.Dopt.y = NaN(1, nbins);
QURU.Dopt.uci = NaN(1, nbins);
QURU.Dopt.lci = NaN(1, nbins);
QURU.S.y = NaN(1, nbins);
QURU.S.uci = NaN(1, nbins);
QURU.S.lci = NaN(1, nbins);

x = quantile(Dstar, 0:(1/nbins):1);
bins = (0+(0.5/nbins)):(1/nbins):(1-(0.5/nbins));
bins = quantile(Dstar, bins);
QURU.Dopt.bins = bins;
QURU.S.bins = bins;

for i = 1:nbins
    
    if i==1
        idx = Dstar>=x(i) & Dstar<=x(i+1);
    else
        idx = Dstar>x(i) & Dstar<=x(i+1);
    end
    Dopt_sub = Dopt_bias(idx);
    S_sub = S_bias(idx);
    
    QURU.Dopt.y(i) = median(Dopt_sub);
    QURU.S.y(i) = median(S_sub);
    
    cis = quantile(Dopt_sub, quants);
    QURU.Dopt.lci(i) = cis(1);
    QURU.Dopt.uci(i) = cis(2);
    
    cis = quantile(S_sub, quants);
    QURU.S.lci(i) = cis(1);
    QURU.S.uci(i) = cis(2);
    
end

% Plots: TSCA
Dstar = tsca(:,1);
Dopt_bias = tsca(:, 6);
S_bias = tsca(:, 9);
S_bias(~isfinite(S_bias)) = NaN;

TSCA.Dopt.y = NaN(1, nbins);
TSCA.Dopt.uci = NaN(1, nbins);
TSCA.Dopt.lci = NaN(1, nbins);
TSCA.S.y = NaN(1, nbins);
TSCA.S.uci = NaN(1, nbins);
TSCA.S.lci = NaN(1, nbins);

x = quantile(Dstar, 0:(1/nbins):1);
bins = (0+(0.5/nbins)):(1/nbins):(1-(0.5/nbins));
bins = quantile(Dstar, bins);
TSCA.Dopt.bins = bins;
TSCA.S.bins = bins;

for i = 1:nbins
    
    if i==1
        idx = Dstar>=x(i) & Dstar<=x(i+1);
    else
        idx = Dstar>x(i) & Dstar<=x(i+1);
    end
    Dopt_sub = Dopt_bias(idx);
    S_sub = S_bias(idx);
    
    TSCA.Dopt.y(i) = median(Dopt_sub);
    TSCA.S.y(i) = median(S_sub);
    
    cis = quantile(Dopt_sub, quants);
    TSCA.Dopt.lci(i) = cis(1);
    TSCA.Dopt.uci(i) = cis(2);
    
    cis = quantile(S_sub, quants);
    TSCA.S.lci(i) = cis(1);
    TSCA.S.uci(i) = cis(2);
    
end



%% Dopt bias
offs = 3;
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 5 5.5];

ax = tight_subplot(5,2,[0.05 0.1],[0.1 0.05], [0.12 0.05]);

axes(ax(1))
fill([PIPO.Dopt.bins fliplr(PIPO.Dopt.bins)], [PIPO.Dopt.lci fliplr(PIPO.Dopt.uci)], clr(1,:), 'FaceAlpha',0.5, 'EdgeColor','none')
set(gca, 'YLim',[-0.44 0.24],'XLim', [0 90], 'XTickLabels','', 'TickDir','out', 'TickLength',[0.025 0.05]);
hold on;
plot([0 90], [0 0], 'k-', 'LineWidth',0.3);
plot(PIPO.Dopt.bins, PIPO.Dopt.y, '-','LineWidth',2, 'Color',clr(1,:).^2);
text(2, 0.19, 'A) \itP. ponderosa', 'FontSize',11);
ylabel('\DeltaD_{opt}^{*} bias (cm)');
box off;

axes(ax(3))
fill([PSME.Dopt.bins fliplr(PSME.Dopt.bins)], [PSME.Dopt.lci fliplr(PSME.Dopt.uci)], clr(2,:), 'FaceAlpha',0.5, 'EdgeColor','none')
set(gca, 'YLim',[-0.44 0.24],'XLim', [0 90], 'XTickLabels','', 'TickDir','out', 'TickLength',[0.025 0.05]);
hold on;
plot([0 90], [0 0], 'k-', 'LineWidth',0.3);
plot(PSME.Dopt.bins, PSME.Dopt.y, '-','LineWidth',2, 'Color',clr(2,:));
text(2, 0.19, 'C) \itP. menziesii', 'FontSize',11);
ylabel('\DeltaD_{opt}^{*} bias (cm)');
box off;

axes(ax(5))
fill([ACRU.Dopt.bins fliplr(ACRU.Dopt.bins)], [ACRU.Dopt.lci fliplr(ACRU.Dopt.uci)], clr(3,:), 'FaceAlpha',0.5, 'EdgeColor','none')
set(gca, 'YLim',[-0.44 0.24],'XLim', [0 90], 'XTickLabels','', 'TickDir','out', 'TickLength',[0.025 0.05]);
hold on;
plot([0 90], [0 0], 'k-', 'LineWidth',0.3);
plot(ACRU.Dopt.bins, ACRU.Dopt.y, '-','LineWidth',2, 'Color',clr(3,:).^2);
text(2, 0.19, 'E) \itA. rubrum', 'FontSize',11);
ylabel('\DeltaD_{opt}^{*} bias (cm)');
box off;

axes(ax(7))
fill([QURU.Dopt.bins fliplr(QURU.Dopt.bins)], [QURU.Dopt.lci fliplr(QURU.Dopt.uci)], clr(4,:), 'FaceAlpha',0.5, 'EdgeColor','none')
set(gca, 'YLim',[-0.44 0.24],'XLim', [0 90], 'XTickLabels','', 'TickDir','out', 'TickLength',[0.025 0.05]);
hold on;
plot([0 90], [0 0], 'k-', 'LineWidth',0.3);
plot(QURU.Dopt.bins, QURU.Dopt.y, '-','LineWidth',2, 'Color',clr(4,:));
text(2, 0.19, 'G) \itQ. rubra', 'FontSize',11);
ylabel('\DeltaD_{opt}^{*} bias (cm)');
box off;

axes(ax(9))
fill([TSCA.Dopt.bins fliplr(TSCA.Dopt.bins)], [TSCA.Dopt.lci fliplr(TSCA.Dopt.uci)], clr(5,:), 'FaceAlpha',0.5, 'EdgeColor','none')
set(gca, 'YLim',[-0.44 0.24],'XLim', [0 90], 'TickDir','out', 'TickLength',[0.025 0.05]);
hold on;
plot([0 90], [0 0], 'k-', 'LineWidth',0.3);
plot(TSCA.Dopt.bins, TSCA.Dopt.y, '-','LineWidth',2, 'Color',clr(5,:).^2);
text(2, 0.19, 'I) \itT. canadensis', 'FontSize',11);
ylabel('\DeltaD_{opt}^{*} bias (cm)');
box off;
xlabel('D^{*} (cm)');


%% S bias
axes(ax(2))
fill([PIPO.S.bins fliplr(PIPO.S.bins)], [PIPO.S.lci fliplr(PIPO.S.uci)], clr(1,:), 'FaceAlpha',0.5, 'EdgeColor','none')
set(gca, 'YLim',[-0.28 1.05],'XLim', [0 90], 'XTickLabels','', 'TickDir','out', 'TickLength',[0.025 0.05]);
hold on;
plot([0 90], [0 0], 'k-', 'LineWidth',0.3);
plot(PIPO.S.bins, PIPO.S.y, '-','LineWidth',2, 'Color',clr(1,:).^2);
text(2+offs, 0.96, 'B) \itP. ponderosa', 'FontSize',11);
ylabel('S^{*} bias');
box off;

axes(ax(4))
fill([PSME.S.bins fliplr(PSME.S.bins)], [PSME.S.lci fliplr(PSME.S.uci)], clr(2,:), 'FaceAlpha',0.5, 'EdgeColor','none')
set(gca, 'YLim',[-0.28 1.05],'XLim', [0 90], 'XTickLabels','', 'TickDir','out', 'TickLength',[0.025 0.05]);
hold on;
plot([0 90], [0 0], 'k-', 'LineWidth',0.3);
plot(PSME.S.bins, PSME.S.y, '-','LineWidth',2, 'Color',clr(2,:));
text(2+offs, 0.96, 'D) \itP. menziesii', 'FontSize',11);
ylabel('S^{*} bias');
box off;

axes(ax(6))
fill([ACRU.S.bins(2:end) fliplr(ACRU.S.bins(2:end))], [ACRU.S.lci(2:end) fliplr(ACRU.S.uci(2:end))], clr(3,:), 'FaceAlpha',0.5, 'EdgeColor','none')
set(gca, 'YLim',[-0.28 1.05],'XLim', [0 90], 'XTickLabels','', 'TickDir','out', 'TickLength',[0.025 0.05]);
hold on;
plot([0 90], [0 0], 'k-', 'LineWidth',0.3);
plot(ACRU.S.bins, ACRU.S.y, '-','LineWidth',2, 'Color',clr(3,:).^2);
text(2+offs, 0.96, 'F) \itA. rubrum', 'FontSize',11);
ylabel('S^{*} bias');
box off;

axes(ax(8))
fill([QURU.S.bins(2:end) fliplr(QURU.S.bins(2:end))], [QURU.S.lci(2:end) fliplr(QURU.S.uci(2:end))], clr(4,:), 'FaceAlpha',0.5, 'EdgeColor','none')
set(gca, 'YLim',[-0.28 1.05],'XLim', [0 90], 'XTickLabels','', 'TickDir','out', 'TickLength',[0.025 0.05]);
hold on;
plot([0 90], [0 0], 'k-', 'LineWidth',0.3);
plot(QURU.S.bins, QURU.S.y, '-','LineWidth',2, 'Color',clr(4,:));
text(2+offs, 0.96, 'H) \itQ. rubra', 'FontSize',11);
ylabel('S^{*} bias');
box off;

axes(ax(10))
fill([TSCA.S.bins(2:end) fliplr(TSCA.S.bins(2:end))], [TSCA.S.lci(2:end) fliplr(TSCA.S.uci(2:end))], clr(5,:), 'FaceAlpha',0.5, 'EdgeColor','none')
set(gca, 'YLim',[-0.28 1.05],'XLim', [0 90], 'TickDir','out', 'TickLength',[0.025 0.05]);
hold on;
plot([0 90], [0 0], 'k-', 'LineWidth',0.3);
plot(TSCA.S.bins, TSCA.S.y, '-','LineWidth',2, 'Color',clr(5,:).^2);
text(2+offs, 0.96, 'J) \itT. canadensis', 'FontSize',11);
ylabel('S^{*} bias');
box off;
xlabel('D^{*} (cm)');

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/dopt_s_bias.tif')
close all;


