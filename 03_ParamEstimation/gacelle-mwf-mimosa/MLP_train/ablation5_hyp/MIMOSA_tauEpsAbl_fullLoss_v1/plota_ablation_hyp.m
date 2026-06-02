%% Plot complex NRMSE versus tau and heatmap ¡ª Nature style
% English comments only

load('ABLA_results_full.mat')

% ---------------- Parse tau and epsA from cfg_name ----------------
nRow     = height(T);
tau_vals = nan(nRow, 1);
eps_vals = nan(nRow, 1);

for i = 1:nRow
    [tau_vals(i), eps_vals(i)] = parse_tau_eps_from_cfg(T.cfg_name(i));
end

T.tau  = tau_vals;
T.epsA = eps_vals;

mask_valid = ~isnan(T.tau) & ~isnan(T.epsA) & ~isnan(T.complex_nrse);
T = T(mask_valid, :);

% ---------------- Convert NRMSE to percentage ----------------
T.complex_nrse_pct = 100 * T.complex_nrse;

tau_unique = unique(T.tau);
eps_unique = unique(T.epsA);

% ---------------- Remove smallest tau ----------------
tau_unique_plot = tau_unique(2:end);   % exclude the minimum tau

% ---------------- Nature line-plot color palette ----------------
nature_colors = [ ...
    0.173, 0.482, 0.714;   % deep steel blue
    0.404, 0.663, 0.812;   % light steel blue
    0.600, 0.780, 0.620;   % muted sage green
    0.650, 0.380, 0.600;   % muted purple   (amber/orange removed)
    0.750, 0.220, 0.220;   % muted red
];

% ---------------- Remove the eps curve that would be orange ----------------
% The 4th curve (index 4) previously got the amber color ¡ª exclude it
% eps_unique_plot = eps_unique([1 2 3 5:end]);   % skip index 4
eps_unique_plot = eps_unique([1 2 3 4 5:end]);   % skip index 4
%% =========================
% Figure 1: line plot
% =========================
fig1 = figure('Color', 'w', 'Position', [100 100 252 200]);
hold on;

legtxt = strings(numel(eps_unique_plot), 1);

for ie = 1:numel(eps_unique_plot)
    eps_i  = eps_unique_plot(ie);
    mask_i = T.epsA == eps_i;
    Ti     = sortrows(T(mask_i, :), 'tau');

    % Exclude the smallest tau row
    Ti = Ti(Ti.tau > min(tau_unique), :);
    Ti = Ti(Ti.tau < max(tau_unique), :);

    c = nature_colors(mod(ie-1, size(nature_colors,1))+1, :);

    plot(Ti.tau, Ti.complex_nrse_pct, '-o', ...
        'Color',           c, ...
        'LineWidth',       2, ...
        'MarkerSize',      8, ...
        'MarkerFaceColor', c, ...
        'MarkerEdgeColor', 'none');

    legtxt(ie) = sprintf('\\epsilon = %.0e', eps_i);
end

% Mark best point (within the filtered data)
T_plot = T(T.tau > min(tau_unique) & ismember(T.epsA, eps_unique_plot), :);
[best_val_pct, best_idx] = min(T_plot.complex_nrse_pct);
best_tau = T_plot.tau(best_idx);
best_eps = T_plot.epsA(best_idx);

plot(best_tau, best_val_pct, 'p', ...
    'MarkerSize',      12, ...
    'LineWidth',       2, ...
    'Color',           [0.750, 0.220, 0.220], ...
    'MarkerFaceColor', [1 1 1]);

% ---------------- Axes ¡ª Nature open-box style ----------------
ax = gca;
ax.FontName      = 'Arial';
ax.FontSize      = 12;
ax.LineWidth     = 1;
ax.TickLength    = [0.02 0.02];
ax.TickDir       = 'out';
ax.Box           = 'off';
ax.Color         = 'none';
ax.XColor        = [0 0 0];
ax.YColor        = [0 0 0];
ax.XGrid         = 'off';
ax.YGrid         = 'off';
ax.GridLineStyle = '-';
ax.GridColor     = [0.85 0.85 0.85];
ax.GridAlpha     = 1.0;
ax.YMinorGrid    = 'off';
ax.XMinorTick    = 'off';
ax.YMinorTick    = 'off';
ax.Layer         = 'top';

xlabel('\tau',              'FontName','Arial','FontSize',12,'Interpreter','tex');
ylabel('Complex NRMSE (%)', 'FontName','Arial','FontSize',12);

lgd = legend(legtxt, ...
    'Location',    'best', ...
    'FontName',    'Arial', ...
    'FontSize',    10, ...
    'Interpreter', 'tex', ...
    'Box',         'off');

set(findall(fig1,'type','text'), 'FontName','Arial','FontSize',12);
hold off;

% ---------------- Export ----------------
outDir = 'fig_tau_eps_nature';
if ~exist(outDir,'dir'), mkdir(outDir); end

exportgraphics(fig1, fullfile(outDir,'Fig_nrmse_vs_tau_nature.png'), 'Resolution', 300);
exportgraphics(fig1, fullfile(outDir,'Fig_nrmse_vs_tau_nature.pdf'), 'ContentType', 'vector');

%% =========================
% Figure 2: heatmap
% =========================
% Build matrix using filtered tau and eps (rows = epsA, cols = tau)
M_pct = nan(numel(eps_unique_plot), numel(tau_unique_plot));

for ie = 1:numel(eps_unique_plot)
    for it = 1:numel(tau_unique_plot)
        mask_ij = (T.epsA == eps_unique_plot(ie)) & (T.tau == tau_unique_plot(it));
        if any(mask_ij)
            M_pct(ie, it) = T.complex_nrse_pct(find(mask_ij, 1, 'first'));
        end
    end
end

fig2 = figure('Color', 'w', 'Position', [100 100 252 180]);

imagesc(M_pct);
axis image;

% White ¡ú deep blue sequential colormap
n_clr = 256;
cmap  = [linspace(0.95, 0.173, n_clr)', ...
         linspace(0.97, 0.482, n_clr)', ...
         linspace(1.00, 0.714, n_clr)'];
colormap(cmap);

cb                = colorbar;
cb.Label.String   = 'Complex NRMSE (%)';
cb.Label.FontName = 'Arial';
cb.Label.FontSize = 12;
cb.FontName       = 'Arial';
cb.FontSize       = 10;
cb.LineWidth      = 0.75;
cb.TickDirection  = 'out';

ax2 = gca;
ax2.FontName  = 'Arial';
ax2.FontSize  = 12;
ax2.LineWidth = 0.75;
ax2.TickDir   = 'out';
ax2.Box       = 'off';
ax2.XColor    = [0 0 0];
ax2.YColor    = [0 0 0];

set(ax2, ...
    'XTick',      1:numel(tau_unique_plot), ...
    'XTickLabel', arrayfun(@(x) sprintf('%.0e',x), tau_unique_plot, 'UniformOutput', false), ...
    'YTick',      1:numel(eps_unique_plot), ...
    'YTickLabel', arrayfun(@(x) sprintf('%.0e',x), eps_unique_plot, 'UniformOutput', false));

xlabel('\tau',     'FontName','Arial','FontSize',12,'Interpreter','tex');
ylabel('\epsilon', 'FontName','Arial','FontSize',12,'Interpreter','tex');

% Cell annotations ¡ª white text on dark, black text on light
vmin = min(M_pct(:), [], 'omitnan');
vmax = max(M_pct(:), [], 'omitnan');
vthr = (vmin + vmax) / 2;

for ie = 1:size(M_pct,1)
    for it = 1:size(M_pct,2)
        if ~isnan(M_pct(ie,it))
            txt_color = 'k';
            if M_pct(ie,it) < vthr
                txt_color = 'w';
            end
            text(it, ie, sprintf('%.2f%%', M_pct(ie,it)), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment',   'middle', ...
                'FontName', 'Arial', ...
                'FontSize', 9, ...
                'Color',    txt_color);
        end
    end
end

% Highlight best cell
[row_best, col_best] = find(M_pct == best_val_pct, 1, 'first');
if ~isempty(row_best)
    rectangle('Position', [col_best-0.5, row_best-0.5, 1, 1], ...
        'EdgeColor', [0 0 0], 'LineWidth', 1.2);
end

set(findall(fig2,'type','text'), 'FontName','Arial','FontSize',10);

exportgraphics(fig2, fullfile(outDir,'Fig_nrmse_heatmap_nature.png'), 'Resolution', 300);
exportgraphics(fig2, fullfile(outDir,'Fig_nrmse_heatmap_nature.pdf'), 'ContentType', 'vector');

%% =========================
% Local function
% =========================
function [tau_val, eps_val] = parse_tau_eps_from_cfg(cfg_name_in)

tau_val = nan;
eps_val = nan;

cfg_str = char(cfg_name_in);
tok = regexp(cfg_str, 'tau_([0-9]+e[m]?[0-9]+)__eps_([0-9]+e[m]?[0-9]+)', 'tokens', 'once');

if isempty(tok), return; end

tau_val = str2double(strrep(tok{1}, 'm', '-'));
eps_val = str2double(strrep(tok{2}, 'm', '-'));

end