clear; clc; close all;

% =========================
% Settings
% =========================
csvFile = 'ABLA_results_sorted.csv';
outDir  = 'fig_temporal_nrmse_bar_FW0_nature';

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
% =========================
T = T(strtrim(T.freq_mode) == "freqRaw", :);
T = T(T.useFW == 0, :);

timeOrder  = ["prepOH", "PT", "PT+prepOH"];
timeLabels = {'Prep one-hot', 'Acquisition index', 'Both'};

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
% Nature-style color palette
% Muted, desaturated: slate blue / steel blue / dusty teal
% =========================
bar_colors = [
    0.404, 0.663, 0.812;   % light steel blue
    0.173, 0.482, 0.714;   % medium blue  (Nature flagship)
    0.071, 0.329, 0.545;   % deep navy blue
];

% =========================
% Figure ¡ª Nature single-column width ¡Ö 89 mm ¡ú ~252 pt
% Use 3.5 in wide; height 2.8 in
% =========================
fig = figure('Color', 'w', 'Position', [100 100 252 200]);   % points ¡Ö 3.5¡Á2.8 in

b = bar(V_nrse_pct, 0.55, 'FaceColor', 'flat');
hold on

for i = 1:numel(V_nrse_pct)
    b.CData(i, :) = bar_colors(i, :);
end

b.EdgeColor = 'none';
b.FaceAlpha = 1.0;

% =========================
% Value labels above bars ¡ª Nature style (small, tight)
% =========================
yRange = max(V_nrse_pct) - min(V_nrse_pct);
yOffset = max(0.04 * yRange, 0.05);

% for i = 1:numel(V_nrse_pct)
%     text(i, V_nrse_pct(i) + yOffset, sprintf('%.2f%%', V_nrse_pct(i)), ...
%         'HorizontalAlignment', 'center', ...
%         'VerticalAlignment', 'bottom', ...
%         'FontSize', 7, ...
%         'FontName', 'Arial', ...
%         'Color', [0.15 0.15 0.15]);
% end

% =========================
% Axes ¡ª Nature style
% Open box (no top/right spine), thin lines, 7 pt font
% =========================
ax = gca;
ax.FontSize   = 7;
ax.FontName   = 'Arial';
ax.LineWidth  = 0.75;
ax.TickLength = [0.02 0.02];
ax.TickDir    = 'out';                  % ticks point outward
ax.Box        = 'off';                  % remove top + right spines
ax.Color      = 'none';                 % transparent axes background
ax.XColor     = [0 0 0];
ax.YColor     = [0 0 0];
ax.XGrid      = 'off';
ax.YGrid      = 'off';
ax.GridLineStyle     = '-';
ax.GridColor         = [0.85 0.85 0.85];
ax.GridAlpha         = 1.0;
ax.YMinorGrid        = 'off';
ax.XMinorTick        = 'off';
ax.YMinorTick        = 'off';
ax.Layer             = 'top';           % grid lines behind bars

set(ax, 'XTick', 1:numel(timeLabels), 'XTickLabel', timeLabels);
xtickangle(0);

xlabel('Temporal feature setting', ...
    'FontSize', 7, 'FontName', 'Arial', 'Color', [0 0 0]);
ylabel('Complex nRMSE (%)', ...
    'FontSize', 7, 'FontName', 'Arial', 'Color', [0 0 0]);

% =========================
% Y-axis limits ¡ª tight, starting at a sensible break
% =========================
ymin = min(V_nrse_pct);
ymax = max(V_nrse_pct);

if ymax > ymin
    padLow  = 0.10 * (ymax - ymin);
    padHigh = 0.25 * (ymax - ymin);
    ylim([max(0, ymin - padLow), ymax + padHigh]);
else
    ylim([0, ymax * 1.3 + 1]);
end

% Force all text objects to Arial
set(findall(gcf, 'type', 'text'), 'FontName', 'Arial', 'FontSize', 12);

hold off

% =========================
% Save ¡ª 300 dpi PNG + vector PDF
% =========================
% exportgraphics(fig, fullfile(outDir, 'Fig_temporal_nRMSE_nature.png'), 'Resolution', 300);
% exportgraphics(fig, fullfile(outDir, 'Fig_temporal_nRMSE_nature.pdf'), 'ContentType', 'vector');
% 
% writetable(summaryTbl, fullfile(outDir, 'temporal_nRMSE_summary_FW0_nature.csv'));