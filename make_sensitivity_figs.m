% New supplemental figures for Global Ecol & Biogeogr paper
clr = wesanderson('fantasticfox1');

%% Bias in DBH measurements (Fig. A12 of diss.)

pipo = csvread('./data/DBH_Comparison_PIPO.csv', 1, 2);
psme = csvread('./data/DBH_Comparison_PSME.csv', 1, 2);
acru = csvread('./data/DBH_Comparison_ACRU.csv', 1, 2);
quco = csvread('./data/DBH_Comparison_QUCO.csv', 1, 2);
quve = csvread('./data/DBH_Comparison_QUVE.csv', 1, 2);
qusp = [quco; quve];

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
plot(acru(:, 2), acru(:, 1), 'o', 'Color',clr(3, :), 'LineWidth',...
    1.2, 'MarkerSize',5);
plot(qusp(:, 2), qusp(:, 1), 'o', 'Color',clr(4, :), 'LineWidth',...
    1.2, 'MarkerSize',5);
[r, p] = corr(pipo, 'type','Spearman');
text(5, 142, '\itP. ponderosa', 'FontSize',11, 'Color',clr(1, :))
text(5, 136, ['n = ', num2str(length(pipo(:, 2)))], 'FontSize',11, 'Color',clr(1, :))
text(5, 130, ['\rho = ', num2str(round(r(1,2), 2))], 'FontSize',11, 'Color',clr(1, :))
[r, p] = corr(psme, 'type','Spearman');
text(5, 118, '\itP. menziesii', 'FontSize',11, 'Color',clr(2, :))
text(5, 112, ['n = ', num2str(length(psme(:, 2)))], 'FontSize',11, 'Color',clr(2, :))
text(5, 106, ['\rho = ', num2str(round(r(1,2), 2))], 'FontSize',11, 'Color',clr(2, :))
[r, p] = corr(acru, 'type','Spearman');
text(45, 142, '\itA. rubrum', 'FontSize',11, 'Color',clr(3, :))
text(45, 136, ['n = ', num2str(length(acru(:, 2)))], 'FontSize',11, 'Color',clr(3, :))
text(45, 130, ['\rho = ', num2str(round(r(1,2), 2))], 'FontSize',11, 'Color',clr(3, :))
[r, p] = corr(qusp, 'type','Spearman');
text(45, 118, '\itQuercus sp.', 'FontSize',11, 'Color',clr(4, :))
text(45, 112, ['n = ', num2str(length(quco(:, 2)))], 'FontSize',11, 'Color',clr(4, :))
text(45, 106, ['\rho = ', num2str(round(r(1,2), 2))], 'FontSize',11, 'Color',clr(4, :))
xlabel('Measured DBH (cm)', 'FontSize',12);
ylabel('Estimated DBH (cm)', 'FontSize',12);
set(gca, 'Position',[0.1300    0.43    0.7750    0.5154]);
text(-23, 150, 'A', 'FontSize',14, 'FontWeight','bold');

subplot(3, 1, 3)
Y1 = histcounts(pipo(:, 1)-pipo(:, 2), -100:5:20);
Y2 = histcounts(psme(:, 1)-psme(:, 2), -100:5:20);
Y3 = histcounts(acru(:, 1)-acru(:, 2), -100:5:20);
Y4 = histcounts(qusp(:, 1)-qusp(:, 2), -100:5:20);
Y = [Y1' Y2' Y3' Y4'];
b = bar((-100+2.5):5:(20-2.5), Y, 1, 'stacked', 'EdgeColor',[0.4 0.4 0.4]);
b(1).FaceColor = clr(1, :);
b(2).FaceColor = clr(2, :);
b(3).FaceColor = clr(3, :);
b(4).FaceColor = clr(4, :);
set(gca, 'YLim',[0 150], 'TickDir','out', 'TickLength',[0.025 0.05], 'FontSize',9);
set(gca, 'Position',[0.1300    0.09    0.7750    0.2157]);
box off;
xlabel('DBH bias (cm)', 'FontSize',12);
ylabel('Count', 'FontSize',12);
text(-124, 150, 'B', 'FontSize',14, 'FontWeight','bold');

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

% QUSP
fprintf('QUSP:\n');
fprintf(['Mean bias: ', num2str(mean(qusp(:, 1)-qusp(:, 2))),'\n']);
fprintf(['Median bias: ', num2str(median(qusp(:, 1)-qusp(:, 2))),'\n']);
fprintf('\n');


%% Sensitivity figures
nbins = 100; % Number of bins within which to calculate median and quantiles
quants = [0.025 0.975];

% 1-D*, 2-D, 3-Dopt*, 4-Dopt, 5-dD, 6-Dopt bias, 7-S*, 8-S, 9-S bias
pipo = csvread('./data/Sensitivity_PIPO.csv', 1, 1);
psme = csvread('./data/Sensitivity_PSME.csv', 1, 1);
acru = csvread('./data/Sensitivity_ACRU.csv', 1, 1);
quco = csvread('./data/Sensitivity_QUCO.csv', 1, 1);
quve = csvread('./data/Sensitivity_QUVE.csv', 1, 1);
qusp = [quco; quve];

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

% Plots: QUSP
Dstar = qusp(:,1);
Dopt_bias = qusp(:, 6);
S_bias = qusp(:, 9);
S_bias(~isfinite(S_bias)) = NaN;

QUSP.Dopt.y = NaN(1, nbins);
QUSP.Dopt.uci = NaN(1, nbins);
QUSP.Dopt.lci = NaN(1, nbins);
QUSP.S.y = NaN(1, nbins);
QUSP.S.uci = NaN(1, nbins);
QUSP.S.lci = NaN(1, nbins);

x = quantile(Dstar, 0:(1/nbins):1);
bins = (0+(0.5/nbins)):(1/nbins):(1-(0.5/nbins));
bins = quantile(Dstar, bins);
QUSP.Dopt.bins = bins;
QUSP.S.bins = bins;

for i = 1:nbins
    
    if i==1
        idx = Dstar>=x(i) & Dstar<=x(i+1);
    else
        idx = Dstar>x(i) & Dstar<=x(i+1);
    end
    Dopt_sub = Dopt_bias(idx);
    S_sub = S_bias(idx);
    
    QUSP.Dopt.y(i) = median(Dopt_sub);
    QUSP.S.y(i) = median(S_sub);
    
    cis = quantile(Dopt_sub, quants);
    QUSP.Dopt.lci(i) = cis(1);
    QUSP.Dopt.uci(i) = cis(2);
    
    cis = quantile(S_sub, quants);
    QUSP.S.lci(i) = cis(1);
    QUSP.S.uci(i) = cis(2);
    
end



%% Dopt bias
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 5 6];

ax = tight_subplot(4,1,0.05,[0.1 0.05], [0.12 0.05]);

axes(ax(1))
fill([PIPO.Dopt.bins fliplr(PIPO.Dopt.bins)], [PIPO.Dopt.lci fliplr(PIPO.Dopt.uci)], clr(1,:), 'FaceAlpha',0.8, 'EdgeColor','none')
set(gca, 'YLim',[-0.44 0.24],'XLim', [0 90], 'XTickLabels','', 'TickDir','out', 'TickLength',[0.015 0.05]);
hold on;
plot([0 90], [0 0], 'k-', 'LineWidth',0.3);
plot(PIPO.Dopt.bins, PIPO.Dopt.y, '-','LineWidth',2, 'Color',clr(1,:)*0.7);
text(2, 0.19, '\itP. ponderosa', 'FontSize',12, 'Color',clr(1,:)*0.7);
ylabel('\DeltaD_{opt}^{*} bias (cm)');
box off;

axes(ax(2))
fill([PSME.Dopt.bins fliplr(PSME.Dopt.bins)], [PSME.Dopt.lci fliplr(PSME.Dopt.uci)], clr(2,:), 'FaceAlpha',0.4, 'EdgeColor','none')
set(gca, 'YLim',[-0.44 0.24],'XLim', [0 90], 'XTickLabels','', 'TickDir','out', 'TickLength',[0.015 0.05]);
hold on;
plot([0 90], [0 0], 'k-', 'LineWidth',0.3);
plot(PSME.Dopt.bins, PSME.Dopt.y, '-','LineWidth',2, 'Color',clr(2,:));
text(2, 0.19, '\itP. menziesii', 'FontSize',12, 'Color',clr(2,:));
ylabel('\DeltaD_{opt}^{*} bias (cm)');
box off;

axes(ax(3))
fill([ACRU.Dopt.bins fliplr(ACRU.Dopt.bins)], [ACRU.Dopt.lci fliplr(ACRU.Dopt.uci)], clr(3,:), 'FaceAlpha',1, 'EdgeColor','none')
set(gca, 'YLim',[-0.44 0.24],'XLim', [0 90], 'XTickLabels','', 'TickDir','out', 'TickLength',[0.015 0.05]);
hold on;
plot([0 90], [0 0], 'k-', 'LineWidth',0.3);
plot(ACRU.Dopt.bins, ACRU.Dopt.y, '-','LineWidth',2, 'Color',clr(3,:)*0.7);
text(2, 0.19, '\itA. rubrum', 'FontSize',12, 'Color',clr(3,:)*0.7);
ylabel('\DeltaD_{opt}^{*} bias (cm)');
box off;

axes(ax(4))
fill([QUSP.Dopt.bins fliplr(QUSP.Dopt.bins)], [QUSP.Dopt.lci fliplr(QUSP.Dopt.uci)], clr(4,:), 'FaceAlpha',0.4, 'EdgeColor','none')
set(gca, 'YLim',[-0.44 0.24],'XLim', [0 90], 'TickDir','out', 'TickLength',[0.015 0.05]);
hold on;
plot([0 90], [0 0], 'k-', 'LineWidth',0.3);
plot(QUSP.Dopt.bins, QUSP.Dopt.y, '-','LineWidth',2, 'Color',clr(4,:));
text(2, 0.19, '\itQuercus sp.', 'FontSize',12, 'Color',clr(4,:));
ylabel('\DeltaD_{opt}^{*} bias (cm)');
box off;
xlabel('D^{*} (cm)');

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/dopt-bias.tif')
close all;


%% S bias
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 5 6];

ax = tight_subplot(4,1,0.05,[0.1 0.05], [0.12 0.05]);

axes(ax(1))
fill([PIPO.S.bins fliplr(PIPO.S.bins)], [PIPO.S.lci fliplr(PIPO.S.uci)], clr(1,:), 'FaceAlpha',0.8, 'EdgeColor','none')
set(gca, 'YLim',[-0.28 1.05],'XLim', [0 90], 'XTickLabels','', 'TickDir','out', 'TickLength',[0.015 0.05]);
hold on;
plot([0 90], [0 0], 'k-', 'LineWidth',0.3);
plot(PIPO.S.bins, PIPO.S.y, '-','LineWidth',2, 'Color',clr(1,:)*0.7);
text(6, 0.96, '\itP. ponderosa', 'FontSize',12, 'Color',clr(1,:)*0.7);
ylabel('S^{*} bias');
box off;

axes(ax(2))
fill([PSME.S.bins fliplr(PSME.S.bins)], [PSME.S.lci fliplr(PSME.S.uci)], clr(2,:), 'FaceAlpha',0.4, 'EdgeColor','none')
set(gca, 'YLim',[-0.28 1.05],'XLim', [0 90], 'XTickLabels','', 'TickDir','out', 'TickLength',[0.015 0.05]);
hold on;
plot([0 90], [0 0], 'k-', 'LineWidth',0.3);
plot(PSME.S.bins, PSME.S.y, '-','LineWidth',2, 'Color',clr(2,:));
text(6, 0.96, '\itP. menziesii', 'FontSize',12, 'Color',clr(2,:));
ylabel('S^{*} bias');
box off;

axes(ax(3))
fill([ACRU.S.bins(2:end) fliplr(ACRU.S.bins(2:end))], [ACRU.S.lci(2:end) fliplr(ACRU.S.uci(2:end))], clr(3,:), 'FaceAlpha',1, 'EdgeColor','none')
set(gca, 'YLim',[-0.28 1.05],'XLim', [0 90], 'XTickLabels','', 'TickDir','out', 'TickLength',[0.015 0.05]);
hold on;
plot([0 90], [0 0], 'k-', 'LineWidth',0.3);
plot(ACRU.S.bins, ACRU.S.y, '-','LineWidth',2, 'Color',clr(3,:)*0.7);
text(6, 0.96, '\itA. rubrum', 'FontSize',12, 'Color',clr(3,:)*0.7);
ylabel('S^{*} bias');
box off;

axes(ax(4))
fill([QUSP.S.bins(2:end) fliplr(QUSP.S.bins(2:end))], [QUSP.S.lci(2:end) fliplr(QUSP.S.uci(2:end))], clr(4,:), 'FaceAlpha',0.4, 'EdgeColor','none')
set(gca, 'YLim',[-0.28 1.05],'XLim', [0 90], 'TickDir','out', 'TickLength',[0.015 0.05]);
hold on;
plot([0 90], [0 0], 'k-', 'LineWidth',0.3);
plot(QUSP.S.bins, QUSP.S.y, '-','LineWidth',2, 'Color',clr(4,:));
text(6, 0.96, '\itQuercus sp.', 'FontSize',12, 'Color',clr(4,:));
ylabel('S^{*} bias');
box off;
xlabel('D^{*} (cm)');

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/s_bias.tif')
close all;


