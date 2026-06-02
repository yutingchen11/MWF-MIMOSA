function T = collect_loss_ablation_results(root_dir, out_csv, out_mat)
% English comments only
%
% Collect loss ablation results directly from saved *_final.mat files.
%
% Example:
%   T = collect_loss_ablation_results( ...
%       'MIMOSA_lossAbl_minChange_v1', ...
%       'ABLA_results_full.csv', ...
%       'ABLA_results_full.mat');

if nargin < 1 || isempty(root_dir)
    root_dir = 'MIMOSA_lossAbl_minChange_v1';
end
if nargin < 2 || isempty(out_csv)
    out_csv = 'ABLA_results_full.csv';
end
if nargin < 3 || isempty(out_mat)
    out_mat = 'ABLA_results_full.mat';
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

        row.cfg_name = get_struct_field_or_default(cfg, 'name', "");
        row.loss_name = get_struct_field_or_default(cfg, 'loss_name', "");
        row.loss_mode = get_struct_field_or_default(cfg, 'loss_mode', NaN);
        row.use_weight = get_struct_field_or_default(cfg, 'use_weight', NaN);
        row.use_norm = get_struct_field_or_default(cfg, 'use_norm', NaN);
        row.tau = get_struct_field_or_default(cfg, 'tau', NaN);
        row.epsA = get_struct_field_or_default(cfg, 'epsA', NaN);

        row.neuron_div = get_struct_field_or_default(cfg, 'neuron_div', NaN);

        numNeurons = get_struct_field_or_default(cfg, 'numNeurons', []);
        row.numNeurons = string(mat2str_safe(numNeurons));

        cfg_feature_idx = get_struct_field_or_default(cfg, 'feature_idx', []);
        row.cfg_feature_idx = string(mat2str_safe(cfg_feature_idx));
    else
        row.cfg_name = "";
        row.loss_name = "";
        row.loss_mode = NaN;
        row.use_weight = NaN;
        row.use_norm = NaN;
        row.tau = NaN;
        row.epsA = NaN;
        row.neuron_div = NaN;
        row.numNeurons = "";
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

        row.mean_abs_gt = get_struct_field_or_default(vm, 'mean_abs_gt', NaN);
        row.mean_abs_pred = get_struct_field_or_default(vm, 'mean_abs_pred', NaN);
        row.mean_abs_err = get_struct_field_or_default(vm, 'mean_abs_err', NaN);

        row.mean_dphi = get_struct_field_or_default(vm, 'mean_dphi', NaN);
        row.std_dphi = get_struct_field_or_default(vm, 'std_dphi', NaN);
        row.mean_abs_dphi = get_struct_field_or_default(vm, 'mean_abs_dphi', NaN);

        row.complex_nrse = get_struct_field_or_default(vm, 'complex_nrse', NaN);
        row.complex_mae = get_struct_field_or_default(vm, 'complex_mae', NaN);

        row.nrse_re = get_struct_field_or_default(vm, 'nrse_re', NaN);
        row.nrse_im = get_struct_field_or_default(vm, 'nrse_im', NaN);

        row.mag_nrse = get_struct_field_or_default(vm, 'mag_nrse', NaN);
        row.mag_mae = get_struct_field_or_default(vm, 'mag_mae', NaN);
    else
        row.mean_abs_gt = NaN;
        row.mean_abs_pred = NaN;
        row.mean_abs_err = NaN;
        row.mean_dphi = NaN;
        row.std_dphi = NaN;
        row.mean_abs_dphi = NaN;
        row.complex_nrse = NaN;
        row.complex_mae = NaN;
        row.nrse_re = NaN;
        row.nrse_im = NaN;
        row.mag_nrse = NaN;
        row.mag_mae = NaN;
    end

    % ---------------- Append ----------------
    rows = [rows; row]; %#ok<AGROW>
end

if isempty(rows)
    warning('No valid *_final.mat files found under: %s', root_dir);
    T = table();
    return;
end

T = struct2table(rows);

% English comment: sort by primary metric
if ismember('complex_nrse', T.Properties.VariableNames)
    T = sortrows(T, 'complex_nrse', 'ascend');
end

% English comment: write outputs
writetable(T, out_csv);
save(out_mat, 'T');

fprintf('Saved CSV: %s\n', out_csv);
fprintf('Saved MAT: %s\n', out_mat);
fprintf('Collected %d runs.\n', height(T));

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