function cfgs = build_ablation_configs_arch_feat()
% English comments only
% Build ablation configs for:
%   - network width: baseNeurons / {1,2,4,8}
%   - feature sets: freqRaw(10-11) OR phaseTrig(12-15) OR both(10-15)
% Constraints:
%   - features 1-9 are mandatory
%   - time features are fixed to PT + prepOH: [16, 17:20]
%   - NO FW feature (22) in any config

% ---- Base neuron setting ----
baseNeurons = [160 240 320 360 480 520 600];

% Divisors for numNeurons ablation
neuronDivs = [1 2 4 8];

% ---- Feature channel mapping (must match your generator) ----
% 1:MWF
% 2:IEW
% 3:T2s MW decay
% 4:T2s IEW decay
% 5:T1 MW norm
% 6:T1 IEW norm
% 7:R2 MW norm
% 8:R2 IEW norm
% 9:B1
% 10:raw Freq MW
% 11:raw Freq IEW
% 12:sin MW
% 13:cos MW
% 14:sin IEW
% 15:cos IEW
% 16:PT
% 17-20:prepOH (4 dims)
% 21:t_acq_norm
% 22:FW (NOT used here)

mandatory = 1:9;

freqA = [10 11];      % raw frequency
freqB = 12:15;        % sin/cos phase trig
freqSets  = {freqA, freqB, [freqA freqB]};
freqNames = {'freqRaw', 'phaseTrig', 'freqRaw+phaseTrig'};

% Fixed time/context set: PT + prepOH
timeSet  = [16, 17:20];
timeName = 'PT+prepOH';

cfgs = struct('id', {}, 'name', {}, 'feature_idx', {}, 'nfeatures', {}, ...
              'freq_mode', {}, 'time_mode', {}, ...
              'neuron_div', {}, 'numNeurons', {});

ic = 0;
for idv = 1:numel(neuronDivs)
    div = neuronDivs(idv);

    % Ensure integer neuron counts
    nn = round(baseNeurons / div);
    nn(nn < 8) = 8;  % safety floor

    for ifr = 1:numel(freqSets)
        ic = ic + 1;

        feat = [mandatory, freqSets{ifr}, timeSet];
        feat = unique(feat, 'stable');

        cfgs(ic).id         = ic; %#ok<AGROW>
        cfgs(ic).freq_mode  = freqNames{ifr}; %#ok<AGROW>
        cfgs(ic).time_mode  = timeName; %#ok<AGROW>
        cfgs(ic).neuron_div = div; %#ok<AGROW>
        cfgs(ic).numNeurons = nn; %#ok<AGROW>
        cfgs(ic).feature_idx = feat; %#ok<AGROW>
        cfgs(ic).nfeatures   = numel(feat); %#ok<AGROW>

        % English comment: stable, filesystem-friendly name
        cfgs(ic).name = sprintf('NDIV%d__%s__%s', div, cfgs(ic).freq_mode, cfgs(ic).time_mode); %#ok<AGROW>
    end
end
end