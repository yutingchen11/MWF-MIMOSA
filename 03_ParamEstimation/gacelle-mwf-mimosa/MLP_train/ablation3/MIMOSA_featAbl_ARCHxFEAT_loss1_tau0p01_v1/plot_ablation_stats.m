% clear; clc; close all;

% =========================
% Settings
% =========================
csvFile = 'ABLA_best_by_group.csv';
outDir  = 'fig_temporal_only_NDIV2';

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

% English comment: convert column names to valid MATLAB identifiers
origNames = T.Properties.VariableNames;
T.Properties.VariableNames = matlab.lang.makeValidName(origNames);

% =========================
% Check required columns
% =========================
requiredVars = {'cfg_name','freq_mode','time_mode','neuron_div', ...
                'complex_nrse','complex_mae','mean_abs_dphi','mag_mae'};

for k = 1:numel(requiredVars)
    assert(ismember(requiredVars{k}, T.Properties.VariableNames), ...
        'Missing required column: %s', requiredVars{k});
end

% =========================
% Convert key columns to string
% =========================
T.cfg_name  = string(T.cfg_name);
T.freq_mode = string(T.freq_mode);
T.time_mode = string(T.time_mode);

% =========================
% Keep freqRaw only
% =========================
freqKey = lower(strtrim(T.freq_mode));
keepFreq = freqKey == "freqraw";
T = T(keepFreq, :);

% =========================
% Keep NDIV = 2 only
% =========================
T = T(T.neuron_div == 2, :);

assert(~isempty(T), 'No rows remain after filtering freqRaw and NDIV=2.');

% =========================
% Normalize temporal labels
% =========================
T.time_label = strings(height(T), 1);
for i = 1:height(T)
    T.time_label(i) = normalizeTimeMode(T.time_mode(i));
end

targetOrder = ["No-time","PT-only","prepOH-only","PT+prepOH"];
keepTime = ismember(T.time_label, targetOrder);
T = T(keepTime, :);

assert(~isempty(T), 'No target temporal modes remain after filtering.');

% =========================
% Keep one row per temporal mode
% English comment: if multiple rows exist, keep the first one
% =========================
Tplot = table();
for j = 1:numel(targetOrder)
    idx = find(T.time_label == targetOrder(j));
    if ~isempty(idx)
        Tplot = [Tplot; T(idx(1), :)]; %#ok<AGROW>
    end
end

assert(~isempty(Tplot), 'No rows available for plotting.');

% =========================
% Compute composite score
% English comment: lower is better for all metrics
% =========================
metrics = {'complex_nrse','complex_mae','mean_abs_dphi','mag_mae'};
scoreMat = zeros(height(Tplot), numel(metrics));

for m = 1:numel(metrics)
    x = Tplot.(metrics{m});
    scoreMat(:, m) = normalizeMinMax(x);
end

Tplot.composite_score = mean(scoreMat, 2);

% =========================
% Print summary
% =========================
disp(' ');
disp('Selected rows for temporal comparison (freqRaw only, NDIV=2):');
disp(Tplot(:, {'time_label','cfg_name','neuron_div','complex_nrse','complex_mae', ...
               'mean_abs_dphi','mag_mae','composite_score'}));

% =========================
% Plot style
% =========================
set(groot, 'defaultAxesFontName', 'Arial');
set(groot, 'defaultTextFontName', 'Arial');
set(groot, 'defaultAxesFontSize', 11);
set(groot, 'defaultAxesLineWidth', 1.0);

timeLabels = cellstr(Tplot.time_label);

% =========================
% Figure 1: grouped bar of four metrics
% =========================
Y = [Tplot.complex_nrse, ...
     Tplot.complex_mae, ...
     Tplot.mean_abs_dphi, ...
     Tplot.mag_mae];

fig1 = figure('Color', 'w', 'Position', [100 100 1100 520]);
bar(Y, 'grouped', 'LineWidth', 1.0);
grid on; box off;

xticks(1:height(Tplot));
xticklabels(timeLabels);
xlabel('Temporal feature setting');
ylabel('Metric value');
title('Effect of temporal features (freqRaw only, NDIV = 2)');
legend({'Complex NRSE','Complex MAE','Mean |¦¤\phi|','Magnitude MAE'}, ...
       'Location', 'northeastoutside', 'Box', 'off');

exportgraphics(fig1, fullfile(outDir, 'Fig1_temporal_grouped_bar_NDIV2.png'), 'Resolution', 300);
exportgraphics(fig1, fullfile(outDir, 'Fig1_temporal_grouped_bar_NDIV2.pdf'), 'ContentType', 'vector');

% =========================
% Figure 2: composite ranking
% =========================
[~, order2] = sort(Tplot.composite_score, 'ascend');
Trank = Tplot(order2, :);

fig2 = figure('Color', 'w', 'Position', [120 120 780 440]);
bar(Trank.composite_score, 'FaceColor', [0.35 0.35 0.35], 'EdgeColor', 'none');
grid on; box off;

xticks(1:height(Trank));
xticklabels(cellstr(Trank.time_label));
ylabel('Composite normalized score');
xlabel('Temporal feature setting');
title('Overall ranking of temporal feature settings (NDIV = 2)');

for i = 1:height(Trank)
    text(i, Trank.composite_score(i), sprintf(' %.3f', Trank.composite_score(i)), ...
        'Rotation', 90, ...
        'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 9);
end

exportgraphics(fig2, fullfile(outDir, 'Fig2_temporal_composite_ranking_NDIV2.png'), 'Resolution', 300);
exportgraphics(fig2, fullfile(outDir, 'Fig2_temporal_composite_ranking_NDIV2.pdf'), 'ContentType', 'vector');

% =========================
% Save summary table
% =========================
writetable(Tplot, fullfile(outDir, 'temporal_summary_NDIV2.csv'));

% =========================
% Local functions
% =========================
function label = normalizeTimeMode(x)
    s = lower(strtrim(string(x)));
    s2 = regexprep(s, '[^a-z0-9]', '');

    if strlength(s2) == 0 || contains(s2, "notime") || contains(s2, "none")
        label = "No-time";
    elseif contains(s2, "pt") && contains(s2, "prepoh")
        label = "PT+prepOH";
    elseif contains(s2, "prepoh") && ~contains(s2, "pt")
        label = "prepOH-only";
    elseif contains(s2, "pt") && ~contains(s2, "prepoh")
        label = "PT-only";
    else
        if contains(s, "pt") && contains(s, "prep")
            label = "PT+prepOH";
        elseif contains(s, "prep")
            label = "prepOH-only";
        elseif contains(s, "pt")
            label = "PT-only";
        else
            label = "No-time";
        end
    end
end

function xNorm = normalizeMinMax(x)
    x = double(x);
    xmin = min(x);
    xmax = max(x);

    if abs(xmax - xmin) < eps
        xNorm = zeros(size(x));
    else
        xNorm = (x - xmin) ./ (xmax - xmin);
    end
end