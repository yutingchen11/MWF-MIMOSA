function T = analyze_ablation_results(root_dir)
% English comments only
% Aggregate feature ablation results and rank by complex_nrse.

if nargin < 1 || isempty(root_dir)
    root_dir = 'MIMOSA_featAbl_freqA_timeOnly_v1';
end

if ~exist(root_dir, 'dir')
    error('Root directory not found: %s', root_dir);
end

out_csv_sorted = fullfile(root_dir, 'ABLA_results_sorted.csv');
out_csv_group  = fullfile(root_dir, 'ABLA_group_best.csv');
out_png        = fullfile(root_dir, 'ABLA_plots.png');

% English comment: list subfolders that look like ablation runs
D = dir(root_dir);
isSub = [D.isdir] & ~ismember({D.name}, {'.','..'});
subDirs = D(isSub);

rows = {};
for i = 1:numel(subDirs)
    run_dir = fullfile(root_dir, subDirs(i).name);

    % English comment: prefer *_final.mat
    f_final = dir(fullfile(run_dir, '*_final.mat'));
    if isempty(f_final)
        f_final = dir(fullfile(run_dir, '*_valMetrics.mat'));
    end

    if isempty(f_final)
        continue;
    end

    fn = fullfile(run_dir, f_final(1).name);
    S = load(fn);

    % English comment: handle both direct and nested valMetrics layouts
    hasVal = isfield(S, 'valMetrics');
    if ~hasVal
        fns = fieldnames(S);
        for k = 1:numel(fns)
            if isstruct(S.(fns{k})) && isfield(S.(fns{k}), 'complex_nrse')
                S.valMetrics = S.(fns{k});
                hasVal = true;
                break;
            end
        end
    end

    if ~hasVal || ~isfield(S.valMetrics, 'complex_nrse')
        continue;
    end

    val = S.valMetrics;

    cfg = struct();
    if isfield(S, 'cfg')
        cfg = S.cfg;
    end

    feature_idx = [];
    if isfield(S, 'feature_idx')
        feature_idx = S.feature_idx;
    end

    final_train_loss = NaN;
    if isfield(S, 'epoch_loss_mean')
        elm = S.epoch_loss_mean;
        if ~isempty(elm)
            final_train_loss = double(elm(end));
        end
    end

    run_name = subDirs(i).name;

    cfg_name  = getfield_default(cfg, 'name', '');
    freq_mode = getfield_default(cfg, 'freq_mode', '');
    time_mode = getfield_default(cfg, 'time_mode', '');

    nfeatures = numel(feature_idx);

    complex_nrse   = getfield_default(val, 'complex_nrse', NaN);
    complex_mae    = getfield_default(val, 'complex_mae', NaN);
    mean_abs_dphi  = getfield_default(val, 'mean_abs_dphi', NaN);
    mag_mae        = getfield_default(val, 'mag_mae', NaN);
    mag_nrse       = getfield_default(val, 'mag_nrse', NaN);
    nrse_re        = getfield_default(val, 'nrse_re', NaN);
    nrse_im        = getfield_default(val, 'nrse_im', NaN);

    rows(end+1,:) = { ...
        run_name, run_dir, fn, cfg_name, freq_mode, time_mode, ...
        nfeatures, mat2str(feature_idx(:).'), final_train_loss, ...
        complex_nrse, complex_mae, mean_abs_dphi, mag_mae, mag_nrse, ...
        nrse_re, nrse_im ...
        }; %#ok<AGROW>
end

if isempty(rows)
    error('No result files found under: %s', root_dir);
end

T = cell2table(rows, 'VariableNames', { ...
    'run_name','run_dir','mat_file','cfg_name','freq_mode','time_mode', ...
    'nfeatures','feature_idx','final_train_loss', ...
    'complex_nrse','complex_mae','mean_abs_dphi','mag_mae','mag_nrse', ...
    'nrse_re','nrse_im'});

% English comment: sort by complex_nrse ascending
T = sortrows(T, 'complex_nrse', 'ascend');

% Save sorted summary
writetable(T, out_csv_sorted);

% Print all results
disp('===== Ranked runs by complex_nrse =====');
disp(T(:, {'run_name','cfg_name','freq_mode','time_mode','nfeatures', ...
           'complex_nrse','complex_mae','mean_abs_dphi','mag_mae'}));

% English comment: keep the best row for each cfg_name
G = findgroups(T.cfg_name);
rowID = (1:height(T)).';

bestRowID = splitapply( ...
    @(nrse, rid) rid(find(nrse == min(nrse), 1, 'first')), ...
    T.complex_nrse, rowID, G);

Tbest = T(bestRowID, :);
Tbest = sortrows(Tbest, 'complex_nrse', 'ascend');

writetable(Tbest, out_csv_group);

% English comment: simple plots
try
    fig = figure('Visible', 'off');
    tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

    % English comment: ranked complex NRSE
    nexttile;
    bar(T.complex_nrse);
    grid on;
    xticks(1:height(T));
    xticklabels(T.time_mode);
    xtickangle(30);
    ylabel('complex\_nrse');
    title('Ranked complex\_nrse');

    % English comment: complex MAE by config
    nexttile;
    bar(T.complex_mae);
    grid on;
    xticks(1:height(T));
    xticklabels(T.time_mode);
    xtickangle(30);
    ylabel('complex\_mae');
    title('complex\_mae by config');

    % English comment: phase error versus complex NRSE
    nexttile;
    scatter(T.complex_nrse, T.mean_abs_dphi, 40, 'filled');
    grid on;
    xlabel('complex\_nrse');
    ylabel('mean|dphi|');
    title('Phase error vs complex\_nrse');

    % English comment: magnitude MAE versus complex NRSE
    nexttile;
    scatter(T.complex_nrse, T.mag_mae, 40, 'filled');
    grid on;
    xlabel('complex\_nrse');
    ylabel('mag\_mae');
    title('Magnitude MAE vs complex\_nrse');

    exportgraphics(fig, out_png);
    close(fig);
catch ME
    warning(ME.identifier, 'Plot export failed: %s', ME.message);
end

fprintf('Saved sorted table: %s\n', out_csv_sorted);
fprintf('Saved group-best  : %s\n', out_csv_group);
fprintf('Saved plots       : %s\n', out_png);

end

function v = getfield_default(S, fn, defaultVal)
% English comments only
if isempty(S) || ~isstruct(S) || ~isfield(S, fn)
    v = defaultVal;
else
    v = S.(fn);
    if isempty(v)
        v = defaultVal;
    end
end
end