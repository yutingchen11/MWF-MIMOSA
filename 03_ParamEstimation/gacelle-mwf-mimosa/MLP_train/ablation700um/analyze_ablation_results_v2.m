function T = analyze_ablation_results_v2(root_dir)
% English comments only
% Robust analyzer for ablation results (supports old/new CSV schemas).

if nargin < 1 || isempty(root_dir)
    error('root_dir is required.');
end
if ~exist(root_dir,'dir')
    error('Root directory not found: %s', root_dir);
end

out_csv_sorted = fullfile(root_dir, 'ABLA_results_sorted_v2.csv');
out_csv_best   = fullfile(root_dir, 'ABLA_best_by_group_v2.csv');

% ---------- Load existing csv if present, otherwise scan folders ----------
csv1 = fullfile(root_dir, 'ABLA_results_sorted.csv');
if exist(csv1,'file')
    T = readtable(csv1);
else
    error('Cannot find %s. Please provide ABLA_results_sorted.csv or run folder scanner.', csv1);
end

% ---------- Helper: safe column getter ----------
hasCol = @(name) any(strcmp(T.Properties.VariableNames, name));

% ---------- Ensure core metric exists ----------
if ~hasCol('complex_nrse')
    error('Table does not contain complex_nrse column.');
end

% ---------- Derive/ensure freq_mode ----------
if ~hasCol('freq_mode')
    % Try parse from cfg_name or run_name
    freq_mode = strings(height(T),1);
    src = "";
    if hasCol('cfg_name'); src = "cfg_name"; end
    if src=="" && hasCol('run_name'); src = "run_name"; end

    if src==""
        freq_mode(:) = "unknown";
    else
        s = string(T.(src));
        % Try patterns used in your naming
        % e.g., "NDIV2__freqRaw__PT+prepOH"
        pat = ["freqRaw\+phaseTrig","freqRaw","phaseTrig"];
        for i = 1:height(T)
            if ~isempty(regexp(s(i), pat(1), 'once'))
                freq_mode(i) = "freqRaw+phaseTrig";
            elseif ~isempty(regexp(s(i), pat(2), 'once'))
                freq_mode(i) = "freqRaw";
            elseif ~isempty(regexp(s(i), pat(3), 'once'))
                freq_mode(i) = "phaseTrig";
            else
                freq_mode(i) = "unknown";
            end
        end
    end
    T.freq_mode = freq_mode;
end

% ---------- Derive/ensure neuron_div ----------
if ~hasCol('neuron_div')
    neuron_div = nan(height(T),1);
    src = "";
    if hasCol('cfg_name'); src="cfg_name"; end
    if src=="" && hasCol('run_name'); src="run_name"; end
    if src==""
        neuron_div(:) = NaN;
    else
        s = string(T.(src));
        for i = 1:height(T)
            tok = regexp(s(i), 'NDIV(\d+)', 'tokens', 'once');
            if ~isempty(tok)
                neuron_div(i) = str2double(tok{1});
            else
                neuron_div(i) = NaN;
            end
        end
    end
    T.neuron_div = neuron_div;
end

% ---------- Ensure time_mode / useFW exist (optional) ----------
if ~hasCol('time_mode')
    % New setup fixes time_mode
    T.time_mode = repmat("PT+prepOH", height(T), 1);
end
if ~hasCol('useFW')
    % New setup has no FW
    T.useFW = zeros(height(T),1);
end

% ---------- Sort by complex_nrse ----------
T = sortrows(T, 'complex_nrse', 'ascend');
writetable(T, out_csv_sorted);

% ---------- Group best: by neuron_div + freq_mode ----------
G = findgroups(T.neuron_div, T.freq_mode);

rowID = (1:height(T))';
bestRowID = splitapply(@(nrse, rid) rid(find(nrse == min(nrse), 1, 'first')), ...
                       T.complex_nrse, rowID, G);

Tbest = T(bestRowID, :);
Tbest = sortrows(Tbest, {'neuron_div','freq_mode'});
writetable(Tbest, out_csv_best);

% ---------- Print summary ----------
disp('===== Best per (neuron_div, freq_mode) =====');
disp(Tbest(:, intersect({'neuron_div','freq_mode','nfeatures','complex_nrse','complex_mae','mean_abs_dphi','mag_mae','run_name','cfg_name'}, Tbest.Properties.VariableNames)));

disp('===== Top-10 overall by complex_nrse =====');
K = min(10,height(T));
disp(T(1:K, intersect({'neuron_div','freq_mode','nfeatures','complex_nrse','complex_mae','mean_abs_dphi','mag_mae','run_name','cfg_name'}, T.Properties.VariableNames)));

fprintf('Saved: %s\n', out_csv_sorted);
fprintf('Saved: %s\n', out_csv_best);
end