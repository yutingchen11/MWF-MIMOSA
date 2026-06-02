clear; clc; close all;

% =========================
% Settings
% =========================
csvFile = 'ABLA_results_sorted.csv';
outDir  = 'fig_feature_nrmse_bar_freqA_brokenAxis';

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
requiredVars = {'freq_mode','time_mode','complex_nrse'};
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
T = T(strtrim(T.freq_mode) == "freqA", :);

timeOrder  = ["none", "prepOH", "PT", "PT+prepOH"];
timeLabels = {'None', 'Prep one-hot', 'Acquisition index', 'Both'};

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
disp('Complex nRMSE (%) for freqA feature ablation:');
summaryTbl = table(timeOrder(:), V_nrse_pct, ...
    'VariableNames', {'time_mode', 'complex_nrmse_percent'});
disp(summaryTbl);

% =========================
% Color palette
% =========================
bar_colors = [
    0.745, 0.827, 0.902;   % light blue
    0.404, 0.663, 0.812;   % steel blue
    0.173, 0.482, 0.714;   % medium blue
    0.071, 0.329, 0.545;   % deep navy
];

% =========================
% Decide whether broken axis is needed
% =========================
valsSorted = sort(V_nrse_pct, 'ascend');

if numel(valsSorted) >= 2
    useBrokenAxis = valsSorted(end) > 4 * valsSorted(end-1);
else
    useBrokenAxis = false;
end

% =========================
% Figure
% =========================
fig = figure('Color', 'w', 'Position', [100 100 320 260]);

if ~useBrokenAxis

    ax = axes('Position', [0.16 0.18 0.78 0.74]);
    b = bar(ax, V_nrse_pct, 0.58, 'FaceColor', 'flat');
    hold(ax, 'on');

    for i = 1:numel(V_nrse_pct)
        b.CData(i, :) = bar_colors(i, :);
    end

    b.EdgeColor = 'none';
    b.FaceAlpha = 1.0;

    ax.FontSize   = 7;
    ax.FontName   = 'Arial';
    ax.LineWidth  = 0.75;
    ax.TickLength = [0.02 0.02];
    ax.TickDir    = 'out';
    ax.Box        = 'off';
    ax.Color      = 'none';
    ax.XColor     = [0 0 0];
    ax.YColor     = [0 0 0];
    ax.XGrid      = 'off';
    ax.YGrid      = 'off';
    ax.XMinorTick = 'off';
    ax.YMinorTick = 'off';

    set(ax, 'XTick', 1:numel(timeLabels), 'XTickLabel', timeLabels);
    xtickangle(ax, 0);

    xlabel(ax, 'Temporal feature setting', ...
        'FontSize', 7, 'FontName', 'Arial', 'Color', [0 0 0]);
    ylabel(ax, 'Complex nRMSE (%)', ...
        'FontSize', 7, 'FontName', 'Arial', 'Color', [0 0 0]);

    ymin = min(V_nrse_pct);
    ymax = max(V_nrse_pct);

    if ymax > ymin
        padLow  = 0.10 * (ymax - ymin);
        padHigh = 0.20 * (ymax - ymin);
        ylim(ax, [max(0, ymin - padLow), ymax + padHigh]);
    else
        ylim(ax, [0, ymax * 1.3 + 1]);
    end

    for i = 1:numel(V_nrse_pct)
        text(ax, i, V_nrse_pct(i) + 0.03 * range(ylim(ax)), ...
            sprintf('%.2f', V_nrse_pct(i)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', 7, ...
            'FontName', 'Arial');
    end

else

    % English comment: split into lower and upper y-ranges
    maxVal = max(V_nrse_pct);
    secondVal = valsSorted(end-1);

    lowMax  = secondVal * 1.35;
    highMin = max(secondVal * 2.0, maxVal * 0.82);
    highMax = maxVal * 1.08;

    % English comment: bottom axis for small values
    ax1 = axes('Position', [0.16 0.16 0.78 0.48]);
    b1 = bar(ax1, V_nrse_pct, 0.58, 'FaceColor', 'flat');
    hold(ax1, 'on');

    for i = 1:numel(V_nrse_pct)
        b1.CData(i, :) = bar_colors(i, :);
    end

    b1.EdgeColor = 'none';
    b1.FaceAlpha = 1.0;

    ax1.FontSize   = 7;
    ax1.FontName   = 'Arial';
    ax1.LineWidth  = 0.75;
    ax1.TickLength = [0.02 0.02];
    ax1.TickDir    = 'out';
    ax1.Box        = 'off';
    ax1.Color      = 'none';
    ax1.XColor     = [0 0 0];
    ax1.YColor     = [0 0 0];
    ax1.XGrid      = 'off';
    ax1.YGrid      = 'off';
    ax1.XMinorTick = 'off';
    ax1.YMinorTick = 'off';

    set(ax1, 'XTick', 1:numel(timeLabels), 'XTickLabel', timeLabels);
    xtickangle(ax1, 0);

    xlabel(ax1, 'Temporal feature setting', ...
        'FontSize', 7, 'FontName', 'Arial', 'Color', [0 0 0]);
    ylabel(ax1, 'Complex nRMSE (%)', ...
        'FontSize', 7, 'FontName', 'Arial', 'Color', [0 0 0]);

    ylim(ax1, [0, lowMax]);

    % English comment: top axis for the very large value
    ax2 = axes('Position', [0.16 0.72 0.78 0.18]);
    b2 = bar(ax2, V_nrse_pct, 0.58, 'FaceColor', 'flat');
    hold(ax2, 'on');

    for i = 1:numel(V_nrse_pct)
        b2.CData(i, :) = bar_colors(i, :);
    end

    b2.EdgeColor = 'none';
    b2.FaceAlpha = 1.0;

    ax2.FontSize   = 7;
    ax2.FontName   = 'Arial';
    ax2.LineWidth  = 0.75;
    ax2.TickLength = [0.02 0.02];
    ax2.TickDir    = 'out';
    ax2.Box        = 'off';
    ax2.Color      = 'none';
    ax2.XColor     = [0 0 0];
    ax2.YColor     = [0 0 0];
    ax2.XGrid      = 'off';
    ax2.YGrid      = 'off';
    ax2.XMinorTick = 'off';
    ax2.YMinorTick = 'off';

    set(ax2, 'XTick', 1:numel(timeLabels), 'XTickLabel', []);
    ylim(ax2, [highMin, highMax]);

    % English comment: add value labels only once, on the upper visible end
    for i = 1:numel(V_nrse_pct)
        if V_nrse_pct(i) <= lowMax
            text(ax1, i, V_nrse_pct(i) + 0.03 * lowMax, ...
                sprintf('%.2f', V_nrse_pct(i)), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', ...
                'FontSize', 7, ...
                'FontName', 'Arial');
        elseif V_nrse_pct(i) >= highMin
            text(ax2, i, V_nrse_pct(i) + 0.03 * (highMax - highMin), ...
                sprintf('%.2f', V_nrse_pct(i)), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', ...
                'FontSize', 7, ...
                'FontName', 'Arial');
        end
    end

    % English comment: draw break marks
    annotation(fig, 'line', [0.145 0.165], [0.655 0.675], 'Color', 'k', 'LineWidth', 0.8);
    annotation(fig, 'line', [0.145 0.165], [0.635 0.655], 'Color', 'k', 'LineWidth', 0.8);
    annotation(fig, 'line', [0.925 0.945], [0.655 0.675], 'Color', 'k', 'LineWidth', 0.8);
    annotation(fig, 'line', [0.925 0.945], [0.635 0.655], 'Color', 'k', 'LineWidth', 0.8);

end

set(findall(gcf, 'type', 'text'), 'FontName', 'Arial');

% =========================
% Save
% =========================
exportgraphics(fig, fullfile(outDir, 'Fig_feature_nRMSE_freqA_brokenAxis.png'), 'Resolution', 300);
exportgraphics(fig, fullfile(outDir, 'Fig_feature_nRMSE_freqA_brokenAxis.pdf'), 'ContentType', 'vector');

writetable(summaryTbl, fullfile(outDir, 'feature_nRMSE_summary_freqA.csv'));