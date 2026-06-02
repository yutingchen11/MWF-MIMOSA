function gen_ablation_list(out_txt, out_csv)
% English comments only
% Generate ablation_list.txt (one cfg name per line) matching run_one_ablation.m.

if nargin < 1 || isempty(out_txt)
    out_txt = 'ablation_list.txt';
end
if nargin < 2 || isempty(out_csv)
    out_csv = 'ablation_list.csv';
end

cfgs = build_feature_ablation_configs();
nCfg = numel(cfgs);

% Write TXT (one per line)
fid = fopen(out_txt, 'w');
if fid < 0
    error('Cannot open %s for writing.', out_txt);
end
for k = 1:nCfg
    fprintf(fid, '%s\n', cfgs(k).name);
end
fclose(fid);

% Write CSV with ID + name
fid = fopen(out_csv, 'w');
if fid < 0
    error('Cannot open %s for writing.', out_csv);
end
fprintf(fid, 'id,name,freq_mode,time_mode,useFW,nfeatures,feature_idx\n');
for k = 1:nCfg
    fi = cfgs(k).feature_idx(:).';
    fprintf(fid, '%d,%s,%s,%s,%d,%d,"%s"\n', ...
        k, cfgs(k).name, cfgs(k).freq_mode, cfgs(k).time_mode, ...
        double(cfgs(k).useFW), numel(fi), strrep(mat2str(fi), '"', ''));
end
fclose(fid);

fprintf('Wrote %s (%d lines)\n', out_txt, nCfg);
fprintf('Wrote %s\n', out_csv);

end

function cfgs = build_feature_ablation_configs()
% English comments only
%
% Full feature channels:
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
% 16:PT index
% 17-20:prepOH (4 dims)
% 21:t_acq_norm
% 22:FW fraction

mandatory = 1:9;

freqA = [10 11];
freqB = 12:15;
freqSets  = {freqA, freqB, [freqA freqB]};
freqNames = {'freqRaw', 'phaseTrig', 'freqRaw+phaseTrig'};

timeGroups = { ...
    [16], ...
    17:20, ...
    [21] ...
};
timeNames = {'PT', 'prepOH', 't_acq'};

% All non-empty subsets of the 3 time groups
timeSets = {};
timeSetNames = {};
for mask = 1:7
    idx = find(bitget(mask, 1:3));
    featTmp = [];
    nameTmp = {};
    for k = 1:numel(idx)
        featTmp = [featTmp timeGroups{idx(k)}]; %#ok<AGROW>
        nameTmp{end+1} = timeNames{idx(k)}; %#ok<AGROW>
    end
    timeSets{end+1} = featTmp; %#ok<AGROW>
    timeSetNames{end+1} = strjoin(nameTmp, '+'); %#ok<AGROW>
end

cfgs = struct('name', {}, 'feature_idx', {}, 'freq_mode', {}, 'time_mode', {}, 'useFW', {});

ic = 0;
for ifr = 1:numel(freqSets)
    for it = 1:numel(timeSets)
        for useFW = 0:1
            ic = ic + 1;

            feat = [mandatory freqSets{ifr} timeSets{it}];
            if useFW == 1
                feat = [feat 22];
            end
            feat = unique(feat, 'stable');

            cfgs(ic).freq_mode = freqNames{ifr}; %#ok<AGROW>
            cfgs(ic).time_mode = timeSetNames{it}; %#ok<AGROW>
            cfgs(ic).useFW     = logical(useFW); %#ok<AGROW>
            cfgs(ic).feature_idx = feat; %#ok<AGROW>
            cfgs(ic).name = sprintf('%s__%s__FW%d', cfgs(ic).freq_mode, cfgs(ic).time_mode, useFW); %#ok<AGROW>
        end
    end
end
end