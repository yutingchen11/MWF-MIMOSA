function [T_runs, T_group, traces] = collect_loss8_epoch_curves(root_dir, out_prefix)
% English comments only
%
% Collect and analyze epoch loss curves from 8-loss ablation runs.
%
% Outputs:
%   T_runs  : per-run table
%   T_group : grouped summary table across seeds
%   traces  : struct containing aligned raw/normalized curves
%
% Example:
%   [T_runs, T_group, traces] = collect_loss8_epoch_curves( ...
%       'MIMOSA_loss8Abl_minChange_v1', ...
%       'loss8_epoch');

if nargin < 1 || isempty(root_dir)
    root_dir = 'MIMOSA_loss8Abl_minChange_v1';
end
if nargin < 2 || isempty(out_prefix)
    out_prefix = 'loss8_epoch';
end

if ~isfolder(root_dir)
    error('Root directory does not exist: %s', root_dir);
end

loss_order = ["WRCL","WCL","RCL","CL","WRL1","WL1","RL1","L1"];

subdirs = dir(root_dir);
subdirs = subdirs([subdirs.isdir]);
subdirs = subdirs(~ismember({subdirs.name}, {'.', '..'}));

rows = [];

for i = 1:numel(subdirs)
    run_dir = fullfile(root_dir, subdirs(i).name);

    final_files = dir(fullfile(run_dir, '*_final.mat'));
    if isempty(final_files)
        fprintf('Skip (no final mat): %s\n', run_dir);
        continue;
    end

    final_file = fullfile(run_dir, final_files(1).name);
    S = load(final_file);

    row = struct();

    % ---------------- Basic identifiers ----------------
    row.run_dir = string(run_dir);
    row.final_file = string(final_file);

    [~, run_name_noext, ~] = fileparts(final_file);
    run_name_noext = erase(run_name_noext, "_final");
    row.run_name = string(run_name_noext);

    % ---------------- Config fields ----------------
    if isfield(S, 'cfg')
        cfg = S.cfg;
        row.cfg_name      = get_struct_field_or_default(cfg, 'name', "");
        row.loss_abbr     = get_struct_field_or_default(cfg, 'loss_abbr', "");
        row.residual_name = get_struct_field_or_default(cfg, 'residual_name', "");
        row.use_weight    = get_struct_field_or_default(cfg, 'use_weight', NaN);
        row.use_norm      = get_struct_field_or_default(cfg, 'use_norm', NaN);
        row.tau           = get_struct_field_or_default(cfg, 'tau', NaN);
        row.epsA          = get_struct_field_or_default(cfg, 'epsA', NaN);
    else
        row.cfg_name      = "";
        row.loss_abbr     = "";
        row.residual_name = "";
        row.use_weight    = NaN;
        row.use_norm      = NaN;
        row.tau           = NaN;
        row.epsA          = NaN;
    end

    % ---------------- Seed ----------------
    if isfield(S, 'seed')
        row.seed = double(S.seed);
    else
        row.seed = NaN;
    end

    % ---------------- Epoch loss curve ----------------
    if isfield(S, 'epoch_loss_mean') && ~isempty(S.epoch_loss_mean)
        loss_curve = double(S.epoch_loss_mean(:));
    else
        fprintf('Skip (no epoch_loss_mean): %s\n', final_file);
        continue;
    end

    row.nepoch = numel(loss_curve);
    row.loss_curve = {loss_curve};

    if isfinite(loss_curve(1)) && abs(loss_curve(1)) > 0
        loss_curve_norm = loss_curve / loss_curve(1);
    else
        loss_curve_norm = nan(size(loss_curve));
    end
    row.loss_curve_norm = {loss_curve_norm};

    % ---------------- Per-run scalar summaries ----------------
    row.initial_loss = loss_curve(1);
    row.final_loss   = loss_curve(end);
    row.best_loss    = min(loss_curve);

    row.initial_loss_norm = loss_curve_norm(1);
    row.final_loss_norm   = loss_curve_norm(end);
    row.best_loss_norm    = min(loss_curve_norm);

    % English comment: normalized total variation
    d1 = diff(loss_curve_norm);
    d2 = diff(loss_curve_norm, 2);

    row.TV_norm = mean(abs(d1), 'omitnan');
    row.Curv_norm = mean(abs(d2), 'omitnan');

    % English comment: fraction of positive steps in normalized curve
    row.IncFrac_norm = mean(d1 > 0, 'omitnan');

    % English comment: epoch to 95% of best improvement
    row.Epoch95 = compute_epoch95(loss_curve_norm);

    % ---------------- Final validation metrics ----------------
    if isfield(S, 'valMetrics')
        vm = S.valMetrics;
        row.complex_nrse  = get_struct_field_or_default(vm, 'complex_nrse', NaN);
        row.complex_mae   = get_struct_field_or_default(vm, 'complex_mae', NaN);
        row.mean_abs_dphi = get_struct_field_or_default(vm, 'mean_abs_dphi', NaN);
        row.mag_mae       = get_struct_field_or_default(vm, 'mag_mae', NaN);
    else
        row.complex_nrse  = NaN;
        row.complex_mae   = NaN;
        row.mean_abs_dphi = NaN;
        row.mag_mae       = NaN;
    end

    rows = [rows; row]; %#ok<AGROW>
end

if isempty(rows)
    warning('No valid runs found under: %s', root_dir);
    T_runs = table();
    T_group = table();
    traces = struct();
    return;
end

T_runs = struct2table(rows);
T_runs.loss_abbr = string(T_runs.loss_abbr);

[~, loss_rank] = ismember(T_runs.loss_abbr, loss_order);
T_runs.loss_rank = loss_rank;
T_runs = sortrows(T_runs, {'loss_rank','seed'}, {'ascend','ascend'});

% ---------------- Save per-run table ----------------
out_csv_runs = sprintf('%s_all_runs.csv', out_prefix);
out_mat_runs = sprintf('%s_all_runs.mat', out_prefix);
writetable(removevars(T_runs, {'loss_curve','loss_curve_norm'}), out_csv_runs);
save(out_mat_runs, 'T_runs');

fprintf('Saved per-run CSV: %s\n', out_csv_runs);
fprintf('Saved per-run MAT: %s\n', out_mat_runs);

% ---------------- Build grouped trace struct ----------------
traces = struct();
group_rows = [];

loss_abbr_unique = loss_order(ismember(loss_order, unique(T_runs.loss_abbr)));

for i = 1:numel(loss_abbr_unique)
    loss_i = loss_abbr_unique(i);
    mask_i = T_runs.loss_abbr == loss_i;
    Ti = T_runs(mask_i, :);

    % English comment: align curves by epoch length
    nRuns = height(Ti);
    maxEpoch = max(Ti.nepoch);

    raw_mat  = nan(nRuns, maxEpoch);
    norm_mat = nan(nRuns, maxEpoch);

    for k = 1:nRuns
        c_raw  = Ti.loss_curve{k};
        c_norm = Ti.loss_curve_norm{k};

        raw_mat(k, 1:numel(c_raw))   = c_raw(:)';
        norm_mat(k, 1:numel(c_norm)) = c_norm(:)';
    end

    traces.(loss_i).raw_mat  = raw_mat;
    traces.(loss_i).norm_mat = norm_mat;
    traces.(loss_i).epoch    = 1:maxEpoch;
    traces.(loss_i).seed     = Ti.seed;

    traces.(loss_i).raw_mean  = mean(raw_mat, 1, 'omitnan');
    traces.(loss_i).raw_std   = std(raw_mat, 0, 1, 'omitnan');
    traces.(loss_i).norm_mean = mean(norm_mat, 1, 'omitnan');
    traces.(loss_i).norm_std  = std(norm_mat, 0, 1, 'omitnan');

    % ---------------- Group summary row ----------------
    grow = struct();
    grow.loss_abbr      = loss_i;
    grow.residual_name  = string(Ti.residual_name(1));
    grow.use_weight     = Ti.use_weight(1);
    grow.use_norm       = Ti.use_norm(1);
    grow.tau            = Ti.tau(1);
    grow.epsA           = Ti.epsA(1);
    grow.n_runs         = nRuns;

    [grow.initial_loss_mean, grow.initial_loss_std] = mean_std_safe(Ti.initial_loss);
    [grow.final_loss_mean,   grow.final_loss_std]   = mean_std_safe(Ti.final_loss);
    [grow.best_loss_mean,    grow.best_loss_std]    = mean_std_safe(Ti.best_loss);

    [grow.final_loss_norm_mean, grow.final_loss_norm_std] = mean_std_safe(Ti.final_loss_norm);
    [grow.best_loss_norm_mean,  grow.best_loss_norm_std]  = mean_std_safe(Ti.best_loss_norm);

    [grow.TV_norm_mean, grow.TV_norm_std]         = mean_std_safe(Ti.TV_norm);
    [grow.Curv_norm_mean, grow.Curv_norm_std]     = mean_std_safe(Ti.Curv_norm);
    [grow.IncFrac_norm_mean, grow.IncFrac_norm_std] = mean_std_safe(Ti.IncFrac_norm);
    [grow.Epoch95_mean, grow.Epoch95_std]         = mean_std_safe(Ti.Epoch95);

    [grow.complex_nrse_mean, grow.complex_nrse_std] = mean_std_safe(Ti.complex_nrse);
    [grow.complex_mae_mean, grow.complex_mae_std]   = mean_std_safe(Ti.complex_mae);
    [grow.mean_abs_dphi_mean, grow.mean_abs_dphi_std] = mean_std_safe(Ti.mean_abs_dphi);
    [grow.mag_mae_mean, grow.mag_mae_std] = mean_std_safe(Ti.mag_mae);

    grow.complex_nrse_mean_std  = format_mean_std(grow.complex_nrse_mean, grow.complex_nrse_std);
    grow.TV_norm_mean_std       = format_mean_std(grow.TV_norm_mean, grow.TV_norm_std);
    grow.Curv_norm_mean_std     = format_mean_std(grow.Curv_norm_mean, grow.Curv_norm_std);
    grow.IncFrac_norm_mean_std  = format_mean_std(grow.IncFrac_norm_mean, grow.IncFrac_norm_std);
    grow.Epoch95_mean_std       = format_mean_std(grow.Epoch95_mean, grow.Epoch95_std);

    group_rows = [group_rows; grow]; %#ok<AGROW>
end

T_group = struct2table(group_rows);

out_csv_group = sprintf('%s_group_summary.csv', out_prefix);
out_mat_group = sprintf('%s_group_summary.mat', out_prefix);
writetable(T_group, out_csv_group);
save(out_mat_group, 'T_group', 'traces');

fprintf('Saved group-summary CSV: %s\n', out_csv_group);
fprintf('Saved group-summary MAT: %s\n', out_mat_group);

% ---------------- Plots ----------------
plot_loss8_mean_std(traces, loss_order, sprintf('%s_norm_loss_curves.png', out_prefix));
plot_pair_mean_std(traces, "WRCL", "WRL1", sprintf('%s_WRCL_vs_WRL1.png', out_prefix));

fprintf('Saved figure: %s\n', sprintf('%s_norm_loss_curves.png', out_prefix));
fprintf('Saved figure: %s\n', sprintf('%s_WRCL_vs_WRL1.png', out_prefix));

end

%% ========================================================================
%% Plot functions
%% ========================================================================

function plot_loss8_mean_std(traces, loss_order, out_png)
% English comments only

fig = figure('Color', 'w', 'Position', [100 100 950 650]);
hold on;

for i = 1:numel(loss_order)
    loss_i = loss_order(i);
    if ~isfield(traces, loss_i)
        continue;
    end

    x = traces.(loss_i).epoch;
    mu = traces.(loss_i).norm_mean;
    sd = traces.(loss_i).norm_std;

    fill([x fliplr(x)], [mu-sd fliplr(mu+sd)], ...
        [0.8 0.8 0.8], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    plot(x, mu, 'LineWidth', 2);
end

ax = gca;
ax.FontName = 'Arial';
ax.FontSize = 13;
ax.FontWeight = 'bold';
ax.LineWidth = 1.6;
ax.Box = 'on';
ax.YGrid = 'on';
ax.XGrid = 'off';
ax.GridAlpha = 0.18;

xlabel('Epoch', 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Normalized train loss', 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold');
title('Mean \pm std of normalized loss curves across seeds', ...
    'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold');

legend(loss_order, 'Location', 'northeastoutside', ...
    'FontName', 'Arial', 'FontSize', 11, 'Box', 'off');

exportgraphics(fig, out_png, 'Resolution', 300);
close(fig);
end

function plot_pair_mean_std(traces, loss_a, loss_b, out_png)
% English comments only

fig = figure('Color', 'w', 'Position', [120 120 900 560]);
hold on;

pair_list = [loss_a, loss_b];
for i = 1:numel(pair_list)
    loss_i = pair_list(i);
    if ~isfield(traces, loss_i)
        continue;
    end

    x = traces.(loss_i).epoch;
    mu = traces.(loss_i).norm_mean;
    sd = traces.(loss_i).norm_std;

    fill([x fliplr(x)], [mu-sd fliplr(mu+sd)], ...
        [0.8 0.8 0.8], 'FaceAlpha', 0.18, 'EdgeColor', 'none');
    plot(x, mu, 'LineWidth', 2.4);
end

ax = gca;
ax.FontName = 'Arial';
ax.FontSize = 13;
ax.FontWeight = 'bold';
ax.LineWidth = 1.6;
ax.Box = 'on';
ax.YGrid = 'on';
ax.XGrid = 'off';
ax.GridAlpha = 0.18;

xlabel('Epoch', 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Normalized train loss', 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold');
title(sprintf('Normalized loss curves: %s vs %s', loss_a, loss_b), ...
    'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold');

legend(pair_list, 'Location', 'best', ...
    'FontName', 'Arial', 'FontSize', 12, 'Box', 'off');

exportgraphics(fig, out_png, 'Resolution', 300);
close(fig);
end

%% ========================================================================
%% Helper functions
%% ========================================================================

function v = get_struct_field_or_default(s, field_name, default_val)
% English comments only

if isstruct(s) && isfield(s, field_name)
    v = s.(field_name);
    if ischar(v)
        v = string(v);
    elseif isempty(v)
        v = default_val;
    end
else
    v = default_val;
end
end

function [m, sd] = mean_std_safe(x)
% English comments only

x = double(x(:));
x = x(isfinite(x));

if isempty(x)
    m = NaN;
    sd = NaN;
else
    m = mean(x);
    sd = std(x, 0);
end
end

function s = format_mean_std(m, sd)
% English comments only

if ~isfinite(m) || ~isfinite(sd)
    s = "";
else
    s = sprintf('%.6g ¡À %.6g', m, sd);
end
end

function e95 = compute_epoch95(loss_curve_norm)
% English comments only

loss_curve_norm = double(loss_curve_norm(:));

if isempty(loss_curve_norm) || any(~isfinite(loss_curve_norm))
    e95 = NaN;
    return;
end

l0 = loss_curve_norm(1);
lbest = min(loss_curve_norm);

if ~isfinite(l0) || ~isfinite(lbest) || abs(l0 - lbest) < eps
    e95 = 1;
    return;
end

target = lbest + 0.05 * (l0 - lbest);
idx = find(loss_curve_norm <= target, 1, 'first');

if isempty(idx)
    e95 = numel(loss_curve_norm);
else
    e95 = idx;
end
end