%% Plot complex NRMSE versus tau and heatmap from ABLA_results_sorted.csv
% English comments only

% clear; clc; close all;

% csv_file = 'ABLA_results_full.csv';
% T = readtable(csv_file, ...
%     'TextType', 'string', ...
%     'VariableNamingRule', 'preserve');
load('ABLA_results_full.mat')

% English comment: parse tau and epsA from cfg_name
nRow = height(T);
tau_vals  = nan(nRow, 1);
eps_vals  = nan(nRow, 1);

for i = 1:nRow
    cfg_name_i = T.cfg_name(i);
    [tau_i, eps_i] = parse_tau_eps_from_cfg(cfg_name_i);
    tau_vals(i) = tau_i;
    eps_vals(i) = eps_i;
end

T.tau  = tau_vals;
T.epsA = eps_vals;

% English comment: keep only rows with valid tau/epsA parsing
mask_valid = ~isnan(T.tau) & ~isnan(T.epsA) & ~isnan(T.complex_nrse);
T = T(mask_valid, :);

% English comment: sort values
tau_unique = unique(T.tau);
eps_unique = unique(T.epsA);

%% =========================
% Figure 1: line plot
% =========================
fig1 = figure('Color', 'w', 'Position', [100 100 820 560]);
hold on;

legtxt = strings(numel(eps_unique), 1);

for ie = 1:numel(eps_unique)
    eps_i = eps_unique(ie);
    mask_i = T.epsA == eps_i;
    Ti = sortrows(T(mask_i, :), 'tau');

    plot(Ti.tau, Ti.complex_nrse, '-o', 'LineWidth', 1.8, 'MarkerSize', 7);
    legtxt(ie) = sprintf('epsA = %.0e', eps_i);
end

% English comment: mark best point
[best_val, best_idx] = min(T.complex_nrse);
best_tau = T.tau(best_idx);
best_eps = T.epsA(best_idx);

plot(best_tau, best_val, 'p', 'MarkerSize', 14, 'LineWidth', 1.8);

text(best_tau, best_val, ...
    sprintf('  best: tau=%.0e, epsA=%.0e, NRMSE=%.6f', best_tau, best_eps, best_val), ...
    'FontSize', 10, 'VerticalAlignment', 'bottom', 'Interpreter', 'none');

set(gca, 'XScale', 'log', 'FontSize', 12, 'LineWidth', 1.0);
xlabel('\tau', 'FontSize', 13);
ylabel('Complex NRMSE', 'FontSize', 13);
title('Complex NRMSE versus \tau under different \epsilon settings', 'FontSize', 14);
legend(legtxt, 'Location', 'best');
grid on;
box on;

exportgraphics(fig1, 'complex_nrmse_vs_tau_by_epsA_matlab.png', 'Resolution', 300);
%%
fig1 = figure('Color', 'w', 'Position', [100 100 820 560]);
hold on;

legtxt = strings(numel(eps_unique), 1);

for ie = 1:numel(eps_unique)
    eps_i = eps_unique(ie);
    mask_i = T.epsA == eps_i;
    Ti = sortrows(T(mask_i, :), 'tau');

    plot(Ti.tau(1:end-1), Ti.complex_nrse(1:end-1), '-o', ...
        'LineWidth', 2.2, ...
        'MarkerSize', 7);
    
    % English comment: use TeX interpreter for epsilon symbol
    legtxt(ie) = sprintf('\\epsilon = %.0e', eps_i);
end

[best_val, best_idx] = min(T.complex_nrse);
best_tau = T.tau(best_idx);
best_eps = T.epsA(best_idx);

plot(best_tau, best_val, 'p', ...
    'MarkerSize', 14, ...
    'LineWidth', 2.0);

% text(best_tau, best_val, ...
%     sprintf('  best: \\tau=%.0e, \\epsilon=%.0e, NRMSE=%.6f', ...
%     best_tau, best_eps, best_val), ...
%     'FontSize', 12, ...
%     'FontName', 'Arial', ...
%     'FontWeight', 'bold', ...
%     'VerticalAlignment', 'bottom', ...
%     'Interpreter', 'tex');

ax = gca;
% ax.XScale = 'log';
ax.FontSize = 14;
ax.FontName = 'Arial';
ax.FontWeight = 'bold';
ax.LineWidth = 1.8;
ax.Box = 'on';

% English comment: make grid visually uniform
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';
ax.GridLineStyle = '-';
ax.GridAlpha = 0.18;
ax.MinorGridAlpha = 0.10;

xlabel('\tau', ...
    'FontSize', 16, ...
    'FontName', 'Arial', ...
    'FontWeight', 'bold', ...
    'Interpreter', 'tex');

ylabel('Complex NRMSE', ...
    'FontSize', 16, ...
    'FontName', 'Arial', ...
    'FontWeight', 'bold');

title('Complex NRMSE versus \tau under different \epsilon settings', ...
    'FontSize', 16, ...
    'FontName', 'Arial', ...
    'FontWeight', 'bold', ...
    'Interpreter', 'tex');

lgd = legend(legtxt, 'Location', 'best');
lgd.FontSize = 12;
lgd.FontName = 'Arial';
lgd.FontWeight = 'bold';
lgd.Interpreter = 'tex';
lgd.Box = 'off';

%% =========================
% Figure 2: heatmap
% =========================
% English comment: build matrix with rows = epsA, cols = tau
M = nan(numel(eps_unique), numel(tau_unique));

for ie = 1:numel(eps_unique)
    for it = 1:numel(tau_unique)
        mask_ij = (T.epsA == eps_unique(ie)) & (T.tau == tau_unique(it));
        if any(mask_ij)
            M(ie, it) = T.complex_nrse(find(mask_ij, 1, 'first'));
        end
    end
end

fig2 = figure('Color', 'w', 'Position', [120 120 820 560]);
imagesc(M);
axis image;

colormap(parula);
cb = colorbar;
cb.Label.String = 'Complex NRMSE';
cb.FontSize = 11;

set(gca, 'XTick', 1:numel(tau_unique), ...
         'XTickLabel', arrayfun(@(x) sprintf('%.0e', x), tau_unique, 'UniformOutput', false), ...
         'YTick', 1:numel(eps_unique), ...
         'YTickLabel', arrayfun(@(x) sprintf('%.0e', x), eps_unique, 'UniformOutput', false), ...
         'FontSize', 12, 'LineWidth', 1.0);

xlabel('\tau', 'FontSize', 13);
ylabel('\epsilon', 'FontSize', 13);
title('Complex NRMSE heatmap over \tau and \epsilon', 'FontSize', 14);
box on;

% English comment: annotate each cell
vmin = min(M(:), [], 'omitnan');
vmax = max(M(:), [], 'omitnan');
vthr = (vmin + vmax) / 2;

for ie = 1:size(M, 1)
    for it = 1:size(M, 2)
        if ~isnan(M(ie, it))
            if M(ie, it) > vthr
                txt_color = 'w';
            else
                txt_color = 'k';
            end

            text(it, ie, sprintf('%.4f', M(ie, it)), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', 10, ...
                'Color', txt_color);
        end
    end
end

% English comment: highlight best cell
[row_best, col_best] = find(M == best_val, 1, 'first');
rectangle('Position', [col_best-0.5, row_best-0.5, 1, 1], ...
          'EdgeColor', 'k', 'LineWidth', 2);

text(col_best, row_best-0.7, 'best', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 10, ...
    'Color', 'k');

exportgraphics(fig2, 'complex_nrmse_heatmap_tau_epsA_matlab.png', 'Resolution', 300);

%% =========================
% Local function
% =========================
function [tau_val, eps_val] = parse_tau_eps_from_cfg(cfg_name_in)
% English comments only

tau_val = nan;
eps_val = nan;

cfg_name_str = char(cfg_name_in);

tok = regexp(cfg_name_str, 'tau_([0-9]+e[m]?[0-9]+)__eps_([0-9]+e[m]?[0-9]+)', 'tokens', 'once');
if isempty(tok)
    return;
end

tau_tag = tok{1};
eps_tag = tok{2};

tau_tag = strrep(tau_tag, 'm', '-');
eps_tag = strrep(eps_tag, 'm', '-');

tau_val = str2double(tau_tag);
eps_val = str2double(eps_tag);

end