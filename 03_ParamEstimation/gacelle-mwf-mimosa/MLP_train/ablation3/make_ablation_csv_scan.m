function T = make_ablation_csv_scan(root_dir)
% English comments only
% Scan all run folders under root_dir and aggregate metrics from *_final.mat.

if nargin < 1 || isempty(root_dir)
    error('root_dir is required.');
end
if ~exist(root_dir,'dir')
    error('Root directory not found: %s', root_dir);
end

out_csv_all  = fullfile(root_dir, 'ABLA_results_scan.csv');
out_csv_best = fullfile(root_dir, 'ABLA_best_by_group.csv');

D = dir(root_dir);
isSub = [D.isdir] & ~ismember({D.name},{'.','..'});
subDirs = D(isSub);

rows = {};
for i = 1:numel(subDirs)
    run_dir = fullfile(root_dir, subDirs(i).name);

    ff = dir(fullfile(run_dir, '*_final.mat'));
    if isempty(ff)
        continue;
    end
    mat_file = fullfile(run_dir, ff(1).name);
    S = load(mat_file);

    if ~isfield(S,'valMetrics') || ~isstruct(S.valMetrics)
        continue;
    end

    vm = S.valMetrics;
    complex_nrse  = getfield_default(vm,'complex_nrse',NaN);
    complex_mae   = getfield_default(vm,'complex_mae',NaN);
    mean_abs_dphi = getfield_default(vm,'mean_abs_dphi',NaN);
    mag_mae       = getfield_default(vm,'mag_mae',NaN);

    final_train_loss = NaN;
    if isfield(S,'epoch_loss_mean') && ~isempty(S.epoch_loss_mean)
        final_train_loss = double(S.epoch_loss_mean(end));
    end

    cfg_name = "";
    freq_mode = "";
    time_mode = "";
    neuron_div = NaN;
    nfeatures = NaN;
    feature_idx_str = "";

    if isfield(S,'cfg') && isstruct(S.cfg)
        cfg_name  = string(getfield_default(S.cfg,'name',""));
        freq_mode = string(getfield_default(S.cfg,'freq_mode',""));
        time_mode = string(getfield_default(S.cfg,'time_mode',""));
        neuron_div = double(getfield_default(S.cfg,'neuron_div',NaN));
    end

    if isfield(S,'feature_idx') && ~isempty(S.feature_idx)
        fi = S.feature_idx(:).';
        nfeatures = numel(fi);
        feature_idx_str = string(mat2str(fi));
    end

    run_name = string(subDirs(i).name);

    % Parse NDIV from run_name if missing
    if isnan(neuron_div)
        tok = regexp(run_name,'NDIV(\d+)','tokens','once');
        if ~isempty(tok), neuron_div = str2double(tok{1}); end
    end

    % Parse freq_mode from run_name if missing
    if strlength(freq_mode)==0 || freq_mode==""
        if contains(run_name,"freqRaw+phaseTrig")
            freq_mode = "freqRaw+phaseTrig";
        elseif contains(run_name,"freqRaw")
            freq_mode = "freqRaw";
        elseif contains(run_name,"phaseTrig")
            freq_mode = "phaseTrig";
        else
            freq_mode = "unknown";
        end
    end

    % Fixed time mode in your new design
    if strlength(time_mode)==0 || time_mode==""
        time_mode = "PT+prepOH";
    end

    rows(end+1,:) = { ...
        run_name, run_dir, mat_file, cfg_name, freq_mode, time_mode, neuron_div, nfeatures, feature_idx_str, ...
        final_train_loss, complex_nrse, complex_mae, mean_abs_dphi, mag_mae ...
        }; %#ok<AGROW>
end

if isempty(rows)
    error('No *_final.mat found under %s (or valMetrics missing).', root_dir);
end

T = cell2table(rows, 'VariableNames', { ...
    'run_name','run_dir','mat_file','cfg_name','freq_mode','time_mode','neuron_div','nfeatures','feature_idx', ...
    'final_train_loss','complex_nrse','complex_mae','mean_abs_dphi','mag_mae'});

% Sort by complex_nrse
T = sortrows(T,'complex_nrse','ascend');
writetable(T, out_csv_all);

% Best per (neuron_div, freq_mode)
G = findgroups(T.neuron_div, T.freq_mode);
rowID = (1:height(T))';
bestRowID = splitapply(@(nrse,rid) rid(find(nrse==min(nrse),1,'first')), T.complex_nrse, rowID, G);
Tbest = T(bestRowID,:);
Tbest = sortrows(Tbest, {'neuron_div','freq_mode','complex_nrse'});
writetable(Tbest, out_csv_best);

disp('===== Top-10 by complex_nrse =====');
K = min(10,height(T));
disp(T(1:K, {'run_name','cfg_name','neuron_div','freq_mode','nfeatures','complex_nrse','complex_mae','mean_abs_dphi','mag_mae'}));

fprintf('Saved: %s\n', out_csv_all);
fprintf('Saved: %s\n', out_csv_best);

end

function v = getfield_default(S, fn, defaultVal)
% English comments only
if isempty(S) || ~isstruct(S) || ~isfield(S, fn)
    v = defaultVal;
else
    v = S.(fn);
    if isempty(v), v = defaultVal; end
end
end