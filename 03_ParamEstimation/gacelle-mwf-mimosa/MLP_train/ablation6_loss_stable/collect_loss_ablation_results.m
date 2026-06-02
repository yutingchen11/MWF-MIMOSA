function [T_runs, T_group] = collect_loss8_ablation_results(root_dir, out_prefix)
% English comments only
%
% Collect 8-loss ablation results directly from saved *_final.mat files.
% This function generates:
%   1) per-run table across all seeds
%   2) grouped summary table (mean/std across seeds for each loss)
%
% Example:
%   [T_runs, T_group] = collect_loss8_ablation_results( ...
%       'MIMOSA_loss8Abl_minChange_v1', ...
%       'loss8_results');

if nargin < 1 || isempty(root_dir)
    root_dir = 'MIMOSA_loss8Abl_minChange_v1';
end
if nargin < 2 || isempty(out_prefix)
    out_prefix = 'loss8_results';
end

if ~isfolder(root_dir)
    error('Root directory does not exist: %s', root_dir);
end

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

        row.cfg_name       = get_struct_field_or_default(cfg, 'name', "");
        row.loss_abbr      = get_struct_field_or_default(cfg, 'loss_abbr', "");
        row.residual_name  = get_struct_field_or_default(cfg, 'residual_name', "");
        row.residual_mode  = get_struct_field_or_default(cfg, 'residual_mode', NaN);
        row.use_weight     = get_struct_field_or_default(cfg, 'use_weight', NaN);
        row.use_norm       = get_struct_field_or_default(cfg, 'use_norm', NaN);
        row.tau            = get_struct_field_or_default(cfg, 'tau', NaN);
        row.epsA           = get_struct_field_or_default(cfg, 'epsA', NaN);
        row.neuron_div     = get_struct_field_or_default(cfg, 'neuron_div', NaN);

        numNeurons = get_struct_field_or_default(cfg, 'numNeurons', []);
        row.numNeurons = string(mat2str_safe(numNeurons));

        cfg_feature_idx = get_struct_field_or_default(cfg, 'feature_idx', []);
        row.cfg_feature_idx = string(mat2str_safe(cfg_feature_idx));
    else
        row.cfg_name       = "";
        row.loss_abbr      = "";
        row.residual_name  = "";
        row.residual_mode  = NaN;
        row.use_weight     = NaN;
        row.use_norm       = NaN;
        row.tau            = NaN;
        row.epsA           = NaN;
        row.neuron_div     = NaN;
        row.numNeurons     = "";
        row.cfg_feature_idx = "";
    end

    % ---------------- Saved feature_idx ----------------
    if isfield(S, 'feature_idx')
        row.feature_idx = string(mat2str_safe(S.feature_idx));
        row.nfeatures = numel(S.feature_idx);
    else
        row.feature_idx = "";
        row.nfeatures = NaN;
    end

    % ---------------- Seed ----------------
    if isfield(S, 'seed')
        row.seed = S.seed;
    else
        row.seed = NaN;
    end

    % ---------------- Final train loss ----------------
    if isfield(S, 'epoch_loss_mean') && ~isempty(S.epoch_loss_mean)
        row.final_train_loss = double(S.epoch_loss_mean(end));
    else
        row.final_train_loss = NaN;
    end

    % ---------------- Validation metrics ----------------
    if isfield(S, 'valMetrics')
        vm = S.valMetrics;

        row.mean_abs_gt   = get_struct_field_or_default(vm, 'mean_abs_gt', NaN);
        row.mean_abs_pred = get_struct_field_or_default(vm, 'mean_abs_pred', NaN);
        row.mean_abs_err  = get_struct_field_or_default(vm, 'mean_abs_err', NaN);

        row.mean_dphi     = get_struct_field_or_default(vm, 'mean_dphi', NaN);
        row.std_dphi      = get_struct_field_or_default(vm, 'std_dphi', NaN);
        row.mean_abs_dphi = get_struct_field_or_default(vm, 'mean_abs_dphi', NaN);

        row.complex_nrse  = get_struct_field_or_default(vm, 'complex_nrse', NaN);
        row.complex_mae   = get_struct_field_or_default(vm, 'complex_mae', NaN);

        row.nrse_re       = get_struct_field_or_default(vm, 'nrse_re', NaN);
        row.nrse_im       = get_struct_field_or_default(vm, 'nrse_im', NaN);

        row.mag_nrse      = get_struct_field_or_default(vm, 'mag_nrse', NaN);
        row.mag_mae       = get_struct_field_or_default(vm, 'mag_mae', NaN);
    else
        row.mean_abs_gt   = NaN;
        row.mean_abs_pred = NaN;
        row.mean_abs_err  = NaN;
        row.mean_dphi     = NaN;
        row.std_dphi      = NaN;
        row.mean_abs_dphi = NaN;
        row.complex_nrse  = NaN;
        row.complex_mae   = NaN;
        row.nrse_re       = NaN;
        row.nrse_im       = NaN;
        row.mag_nrse      = NaN;
        row.mag_mae       = NaN;
    end

    rows = [rows; row]; %#ok<AGROW>
end

if isempty(rows)
    warning('No valid *_final.mat files found under: %s', root_dir);
    T_runs = table();
    T_group = table();
    return;
end

T_runs = struct2table(rows);

% English comment: custom order of loss abbreviations
loss_order = ["WRCL","WCL","RCL","CL","WRL1","WL1","RL1","L1"];
T_runs.loss_abbr = string(T_runs.loss_abbr);

% English comment: sort per-run table
[~, loss_rank] = ismember(T_runs.loss_abbr, loss_order);
T_runs.loss_rank = loss_rank;
T_runs = sortrows(T_runs, {'loss_rank','seed','complex_nrse'}, {'ascend','ascend','ascend'});

% English comment: save per-run table
out_csv_runs = sprintf('%s_all_runs.csv', out_prefix);
out_mat_runs = sprintf('%s_all_runs.mat', out_prefix);

writetable(T_runs, out_csv_runs);
save(out_mat_runs, 'T_runs');

fprintf('Saved per-run CSV: %s\n', out_csv_runs);
fprintf('Saved per-run MAT: %s\n', out_mat_runs);

% ---------------- Group summary across seeds ----------------
loss_abbr_unique = loss_order(ismember(loss_order, unique(T_runs.loss_abbr)));

group_rows = [];

for i = 1:numel(loss_abbr_unique)
    loss_i = loss_abbr_unique(i);
    mask_i = T_runs.loss_abbr == loss_i;
    Ti = T_runs(mask_i, :);

    grow = struct();
    grow.loss_abbr = loss_i;

    if any(mask_i)
        grow.residual_name = string(Ti.residual_name(1));
        grow.use_weight    = Ti.use_weight(1);
        grow.use_norm      = Ti.use_norm(1);
        grow.tau           = Ti.tau(1);
        grow.epsA          = Ti.epsA(1);
        grow.nfeatures     = Ti.nfeatures(1);
        grow.feature_idx   = string(Ti.feature_idx(1));
        grow.n_runs        = height(Ti);

        % English comment: compute mean/std across seeds
        [grow.final_train_loss_mean, grow.final_train_loss_std] = mean_std_safe(Ti.final_train_loss);

        [grow.complex_nrse_mean, grow.complex_nrse_std] = mean_std_safe(Ti.complex_nrse);
        [grow.complex_mae_mean,  grow.complex_mae_std]  = mean_std_safe(Ti.complex_mae);
        [grow.mean_abs_dphi_mean, grow.mean_abs_dphi_std] = mean_std_safe(Ti.mean_abs_dphi);
        [grow.mag_mae_mean, grow.mag_mae_std] = mean_std_safe(Ti.mag_mae);

        [grow.nrse_re_mean, grow.nrse_re_std] = mean_std_safe(Ti.nrse_re);
        [grow.nrse_im_mean, grow.nrse_im_std] = mean_std_safe(Ti.nrse_im);
        [grow.mag_nrse_mean, grow.mag_nrse_std] = mean_std_safe(Ti.mag_nrse);

        % English comment: store printable mean¡Àstd strings
        grow.complex_nrse_mean_std = format_mean_std(grow.complex_nrse_mean, grow.complex_nrse_std);
        grow.complex_mae_mean_std = format_mean_std(grow.complex_mae_mean, grow.complex_mae_std);
        grow.mean_abs_dphi_mean_std = format_mean_std(grow.mean_abs_dphi_mean, grow.mean_abs_dphi_std);
        grow.mag_mae_mean_std = format_mean_std(grow.mag_mae_mean, grow.mag_mae_std);
    else
        grow.residual_name = "";
        grow.use_weight = NaN;
        grow.use_norm = NaN;
        grow.tau = NaN;
        grow.epsA = NaN;
        grow.nfeatures = NaN;
        grow.feature_idx = "";
        grow.n_runs = 0;

        grow.final_train_loss_mean = NaN;
        grow.final_train_loss_std = NaN;

        grow.complex_nrse_mean = NaN;
        grow.complex_nrse_std = NaN;
        grow.complex_mae_mean = NaN;
        grow.complex_mae_std = NaN;
        grow.mean_abs_dphi_mean = NaN;
        grow.mean_abs_dphi_std = NaN;
        grow.mag_mae_mean = NaN;
        grow.mag_mae_std = NaN;

        grow.nrse_re_mean = NaN;
        grow.nrse_re_std = NaN;
        grow.nrse_im_mean = NaN;
        grow.nrse_im_std = NaN;
        grow.mag_nrse_mean = NaN;
        grow.mag_nrse_std = NaN;

        grow.complex_nrse_mean_std = "";
        grow.complex_mae_mean_std = "";
        grow.mean_abs_dphi_mean_std = "";
        grow.mag_mae_mean_std = "";
    end

    group_rows = [group_rows; grow]; %#ok<AGROW>
end

T_group = struct2table(group_rows);

out_csv_group = sprintf('%s_group_summary.csv', out_prefix);
out_mat_group = sprintf('%s_group_summary.mat', out_prefix);

writetable(T_group, out_csv_group);
save(out_mat_group, 'T_group');

fprintf('Saved group-summary CSV: %s\n', out_csv_group);
fprintf('Saved group-summary MAT: %s\n', out_mat_group);
fprintf('Collected %d runs across %d loss types.\n', height(T_runs), height(T_group));

end

%% ========================================================================
%% Local helpers
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

function s = mat2str_safe(x)
% English comments only

if isempty(x)
    s = '';
    return;
end

if isstring(x)
    s = char(x);
    return;
end

if ischar(x)
    s = x;
    return;
end

try
    s = mat2str(x);
catch
    s = '';
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