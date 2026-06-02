clear; clc; close all;

% =========================
% Settings
% =========================
csvFile = 'ABLA_results_sorted.csv';
outDir  = 'fig_temporal_nrmse_bar_FW0_pretty';

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% =========================
% Read table
% =========================
opts = detectImportOptions(csvFile, ...
    'Delimiter', ',', ...
    'VariableNamingRule', 'preserve');

T = readtable(csvFile, opts);

% English comment: convert original column names to valid MATLAB identifiers
T.Properties.VariableNames = matlab.lang.makeValidName(T.Properties.VariableNames);

% =========================
% Check required columns
% =========================
requiredVars = {'freq_mode','time_mode','useFW','complex_nrse'};
for k = 1:numel(requiredVars)
    assert(ismember(requiredVars{k}, T.Properties.VariableNames), ...
        'Missing required column: %s', requiredVars{k});
end

% =========================
% Convert key columns to string
% =========================
T.freq_mode = string(T.freq_mode);
T.time_mode = string(T.time_mode);

% =========================
% Filter data
% 1) freqRaw only
% 2) FW = 0 only
% 3) keep only prepOH / PT / PT+prepOH
% =========================
T = T(strtrim(T.freq_mode) == "freqRaw", :);
T = T(T.useFW == 0, :);

timeOrder  = ["prepOH", "PT", "PT+prepOH"];
timeLabels = {'prep one-hot', 'acquisition index', 'both'};

T = T(ismember(strtrim(T.time_mode), timeOrder), :);
assert(~isempty(T), 'No rows remain after filtering.');

% =========================
% Build nRMSE vector
% =========================
V_nrse = nan(numel(timeOrder), 1);

for i = 1:numel(timeOrder)
    idx = strtrim(T.time_mode) == timeOrder(i);
    assert(any(idx), 'Missing time_mode = %s', timeOrder(i));
    vals = T.complex_nrse(idx);
    V_nrse(i) = vals(1);
end

% English comment: convert to percentage
V_nrse_pct = 100 * V_nrse(:);

% =========================
% Print summary
% =========================
disp(' ');
disp('Complex nRMSE (%) for freqRaw only, FW = 0:');
summaryTbl = table(timeOrder(:), V_nrse_pct, ...
    'VariableNames', {'time_mode', 'complex_nrmse_percent'});
disp(summaryTbl);

% =========================
% Plot style
% =========================
base_color = [0.961 0.651 0.137];   % English comment: amber

% English comment: create three shades from the base color
bar_colors = [
    min(base_color + 0.18, 1.0);
    base_color;
    max(base_color - 0.12, 0.0)
];

fig = figure('Color', 'w', 'Position', [100 100 620 430]);

b = bar(V_nrse_pct, 0.62, 'FaceColor', 'flat');
hold on

for i = 1:numel(V_nrse_pct)
    b.CData(i, :) = bar_colors(i, :);
end

b.EdgeColor = 'none';
b.FaceAlpha = 0.90;

% % English comment: highlight the best setting
% [bestVal, bestIdx] = min(V_nrse_pct);
% plot(bestIdx, bestVal, 'ko', 'MarkerSize', 8, 'LineWidth', 1.2);

% English comment: add value labels above bars
yRange = max(V_nrse_pct) - min(V_nrse_pct);
if yRange < 1e-8
    yOffset = 0.3;
else
    yOffset = 0.08 * yRange;
end

for i = 1:numel(V_nrse_pct)
    text(i, V_nrse_pct(i) + yOffset, sprintf('%.2f%%', V_nrse_pct(i)), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 11, ...
        'FontName', 'Arial');
end

% =========================
% Axes formatting
% =========================
set(gca, 'XTick', 1:numel(timeLabels), 'XTickLabel', timeLabels);

xlabel('Temporal feature setting', ...
    'FontSize', 14, 'FontName', 'Arial');

ylabel('Complex nRMSE (%)', ...
    'FontSize', 14,  'FontName', 'Arial');

% title('Effect of temporal features on complex nRMSE', ...
%     'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');

ax = gca;
ax.FontSize = 14;
ax.FontName = 'Arial';
ax.LineWidth = 1.5;
ax.Box = 'off';
ax.Color = 'w';
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.XGrid = 'off';
ax.YGrid = 'off';

xtickangle(20);

% English comment: make the y-axis tighter for better comparison
ymin = min(V_nrse_pct);
ymax = max(V_nrse_pct);

if ymax > ymin
    padLow  = 0.12 * (ymax - ymin);
    padHigh = 0.22 * (ymax - ymin);
    ylim([max(0, ymin - padLow), ymax + padHigh]);
else
    ylim([0, ymax * 1.3 + 1]);
end

set(findall(gcf, 'type', 'text'), 'FontName', 'Arial');

hold off

% =========================
% Save outputs
% =========================
% exportgraphics(fig, fullfile(outDir, 'Fig_temporal_complex_nRMSE_bar_FW0_pretty.png'), 'Resolution', 300);
% exportgraphics(fig, fullfile(outDir, 'Fig_temporal_complex_nRMSE_bar_FW0_pretty.pdf'), 'ContentType', 'vector');
% 
% writetable(summaryTbl, fullfile(outDir, 'temporal_complex_nRMSE_summary_FW0_pretty.csv'));