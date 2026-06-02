function T = analyze_ablation_results(root_dir)
% English comments only
% Aggregate ablation results and rank by complex_nrse.

if nargin < 1 || isempty(root_dir)
    root_dir = 'MIMOSA_tauEpsAbl_fullLoss_v1';
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

    % English comment: prefer *_final.mat produced by run_one_ablation
    f_final = dir(fullfile(run_dir, '*_final.mat'));
    if isempty(f_final)
        % fallback: maybe *_valMetrics.mat
        f_final = dir(fullfile(run_dir, '*_valMetrics.mat'));
    end

    if isempty(f_final)
        % English comment: skip incomplete runs
        continue;
    end

    fn = fullfile(run_dir, f_final(1).name);

    S = load(fn);

    % English comment: handle both *_final.mat and *_valMetrics.mat layouts
    hasVal = isfield(S, 'valMetrics');
    if ~hasVal
        % try nested fields (unlikely, but safe)
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

    % English comment: extract common metadata
    val = S.valMetrics;

    cfg = struct();
    if isfield(S, 'cfg'); cfg = S.cfg; end

    feature_idx = [];
    if isfield(S, 'feature_idx'); feature_idx = S.feature_idx; end

    final_train_loss = NaN;
    if isfield(S, 'epoch_loss_mean')
        elm = S.epoch_loss_mean;
        if ~isempty(elm)
            final_train_loss = double(elm(end));
        end
    end

    % English comment: parse run name from folder
    run_name = subDirs(i).name;

    % English comment: safe getters
    cfg_name  = getfield_default(cfg, 'name', '');
    freq_mode = getfield_default(cfg, 'freq_mode', '');
    time_mode = getfield_default(cfg, 'time_mode', '');
    useFW     = getfield_default(cfg, 'useFW', NaN);

    if islogical(useFW)
        useFW = double(useFW);
    end

    nfeatures = numel(feature_idx);

    % English comment: store metrics (main focus: complex_nrse)
    complex_nrse = getfield_default(val, 'complex_nrse', NaN);
    complex_mae  = getfield_default(val, 'complex_mae',  NaN);
    mean_abs_dphi = getfield_default(val, 'mean_abs_dphi', NaN);
    mag_mae      = getfield_default(val, 'mag_mae', NaN);

    rows(end+1,:) = { ...
        run_name, run_dir, fn, cfg_name, freq_mode, time_mode, useFW, nfeatures, mat2str(feature_idx(:).'), ...
        final_train_loss, complex_nrse, complex_mae, mean_abs_dphi, mag_mae ...
        }; %#ok<AGROW>
end

if isempty(rows)
    error('No result files found under: %s', root_dir);
end

T = cell2table(rows, 'VariableNames', { ...
    'run_name','run_dir','mat_file','cfg_name','freq_mode','time_mode','useFW','nfeatures','feature_idx', ...
    'final_train_loss','complex_nrse','complex_mae','mean_abs_dphi','mag_mae'});

% English comment: sort by complex_nrse (ascending)
T = sortrows(T, 'complex_nrse', 'ascend');

% Save sorted summary
writetable(T, out_csv_sorted);

% Print top-10
K = min(10, height(T));
disp('===== Top runs by complex_nrse =====');
disp(T(1:K, {'run_name','cfg_name','freq_mode','time_mode','useFW','nfeatures','complex_nrse','complex_mae','mean_abs_dphi','mag_mae'}));

% English comment: group-wise best (freq_mode + time_mode + useFW)
% English comments only
T = readtable(fullfile(root_dir,'ABLA_results_sorted.csv'));

% Group by keys you want
G = findgroups(T.freq_mode, T.time_mode, T.useFW);

rowID = (1:height(T))';  % global row index

bestRowID = splitapply(@(nrse, rid) rid(find(nrse == min(nrse), 1, 'first')), ...
                       T.complex_nrse, rowID, G);

Tbest = T(bestRowID, :);
Tbest = sortrows(Tbest, 'complex_nrse', 'ascend');

writetable(Tbest, fullfile(root_dir, 'ABLA_group_best_fixed.csv'));

% Simple plots
try
    fig = figure('Visible','off');
    tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

    nexttile;
    plot(T.complex_nrse, '-o'); grid on;
    xlabel('Rank'); ylabel('complex\_nrse');
    title('Ranked complex\_nrse');

    nexttile;
    scatter(T.nfeatures, T.complex_nrse, 20, 'filled'); grid on;
    xlabel('nfeatures'); ylabel('complex\_nrse');
    title('complex\_nrse vs nfeatures');

    nexttile;
    scatter(T.complex_nrse, T.mean_abs_dphi, 20, 'filled'); grid on;
    xlabel('complex\_nrse'); ylabel('mean|dphi|');
    title('Phase error vs complex\_nrse');

    nexttile;
    scatter(T.complex_nrse, T.mag_mae, 20, 'filled'); grid on;
    xlabel('complex\_nrse'); ylabel('mag\_mae');
    title('Magnitude MAE vs complex\_nrse');

    exportgraphics(fig, out_png);
    close(fig);
catch ME
    warning(E.identifier,'Plot export failed: %s', ME.message);
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