%% Plot complex NRMSE for different loss settings ¡ª Nature style
clear; clc; close all;

% ---------------- Data ----------------
cfg_name = { ...
    'LOSS01__relCplx__W1__N1', ...
    'LOSS02__relCplx__W1__N0', ...
    'LOSS05__plainL1', ...
    'LOSS04__relCplx__W0__N0', ...
    'LOSS03__relCplx__W0__N1', ...
    'LOSS06__plainL2'};

complex_nrse = [ ...
    0.0016754452, ...
    0.0021368645, ...
    0.0022423828, ...
    0.0024065680, ...
    0.0032924288, ...
    0.0078818863 ];

xticklabels_short = {'WRCL','WCL','L1','CL','RCL','L2'};

% ---------------- Sort ascending ----------------
[complex_nrse_sorted, idx] = sort(complex_nrse, 'ascend');
cfg_name_sorted      = cfg_name(idx);
xticklabels_sorted   = xticklabels_short(idx);

% ---------------- Nature-palette category colors ----------------
% Muted, desaturated ¡ª consistent with Nature figure guidelines
color_full   = [0.173, 0.482, 0.714];   % deep steel blue   (best / full loss)
color_relcpx = [0.404, 0.663, 0.812];   % light steel blue  (other relCplx)
color_plain  = [0.600, 0.780, 0.620];   % muted sage green  (plain L1/L2)

bar_colors = zeros(numel(cfg_name_sorted), 3);
for i = 1:numel(cfg_name_sorted)
    name_i = cfg_name_sorted{i};
    if strcmp(name_i, 'LOSS01__relCplx__W1__N1')
        bar_colors(i,:) = color_full;
    elseif contains(name_i, 'relCplx')
        bar_colors(i,:) = color_relcpx;
    elseif contains(name_i, 'plainL1') || contains(name_i, 'plainL2')
        bar_colors(i,:) = color_plain;
    else
        bar_colors(i,:) = [0.65 0.65 0.65];
    end
end

% ---------------- Figure ¡ª Nature single-column width ¡Ö 89 mm ----------------
fig = figure('Color', 'w', 'Position', [100 100 252 200]);

b = bar(complex_nrse_sorted, 0.55, 'FaceColor', 'flat');
hold on;

for i = 1:numel(complex_nrse_sorted)
    b.CData(i,:) = bar_colors(i,:);
end
b.EdgeColor = 'none';
b.FaceAlpha = 1.0;

% ---------------- Value labels ----------------
yRange  = max(complex_nrse_sorted) - min(complex_nrse_sorted);
yOffset = max(0.03 * yRange, 5e-5);

for i = 1:numel(complex_nrse_sorted)
    text(i, complex_nrse_sorted(i) + yOffset, ...
        sprintf('%.4f', complex_nrse_sorted(i)), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment',   'bottom', ...
        'FontName', 'Arial', ...
        'FontSize', 5.5, ...
        'Color',    [0.15 0.15 0.15]);
end

% ---------------- Axes ¡ª Nature open-box style ----------------
ax = gca;
ax.FontName          = 'Arial';
ax.FontSize          = 7;
ax.LineWidth         = 0.75;
ax.TickLength        = [0.02 0.02];
ax.TickDir           = 'out';
ax.Box               = 'off';
ax.Color             = 'none';
ax.XColor            = [0 0 0];
ax.YColor            = [0 0 0];
ax.XGrid             = 'off';
ax.YGrid             = 'on';
ax.GridLineStyle     = '-';
ax.GridColor         = [0.85 0.85 0.85];
ax.GridAlpha         = 1.0;
ax.YMinorGrid        = 'off';
ax.XMinorTick        = 'off';
ax.YMinorTick        = 'off';
ax.Layer             = 'top';

ax.XTick             = 1:numel(xticklabels_sorted);
ax.XTickLabel        = xticklabels_sorted;
ax.XTickLabelRotation = 0;

xlabel('Loss setting',   'FontName','Arial','FontSize',7,'Color',[0 0 0]);
ylabel('Complex NRMSE', 'FontName','Arial','FontSize',7,'Color',[0 0 0]);

% Tight y-limits
ymin = min(complex_nrse_sorted);
ymax = max(complex_nrse_sorted);
padLow  = 0.10 * (ymax - ymin);
padHigh = 0.28 * (ymax - ymin);
ylim([max(0, ymin - padLow), ymax + padHigh]);

% ---------------- Legend ----------------
h1 = bar(nan, 'FaceColor', color_full,   'EdgeColor', 'none');
h2 = bar(nan, 'FaceColor', color_relcpx, 'EdgeColor', 'none');
h3 = bar(nan, 'FaceColor', color_plain,  'EdgeColor', 'none');

lgd = legend([h1 h2 h3], ...
    {'Full loss', 'Other relCplx', 'Plain L1/L2'}, ...
    'Location',  'northeast', ...
    'FontName',  'Arial', ...
    'FontSize',  6, ...
    'Box',       'off');

% Force all text to Arial 7 pt
set(findall(gcf,'type','text'), 'FontName','Arial','FontSize',7);

hold off;

% ---------------- Export ----------------
outDir = 'fig_loss_nrmse_nature';
if ~exist(outDir,'dir'), mkdir(outDir); end

exportgraphics(fig, fullfile(outDir,'Fig_loss_nrmse_nature.png'), 'Resolution', 300);
exportgraphics(fig, fullfile(outDir,'Fig_loss_nrmse_nature.pdf'), 'ContentType', 'vector');