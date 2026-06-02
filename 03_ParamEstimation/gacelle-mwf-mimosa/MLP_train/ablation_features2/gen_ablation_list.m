function gen_ablation_list(out_txt, out_csv)
% English comments only
% Generate ablation_list.txt (one cfg name per line).

if nargin < 1 || isempty(out_txt)
    out_txt = 'ablation_list.txt';
end
if nargin < 2 || isempty(out_csv)
    out_csv = 'ablation_list.csv';
end

cfgs = build_feature_ablation_configs();
nCfg = numel(cfgs);

% Write TXT (one cfg name per line)
fid = fopen(out_txt, 'w');
if fid < 0
    error('Cannot open %s for writing.', out_txt);
end
for k = 1:nCfg
    fprintf(fid, '%s\n', cfgs(k).name);
end
fclose(fid);

% Write CSV with metadata
fid = fopen(out_csv, 'w');
if fid < 0
    error('Cannot open %s for writing.', out_csv);
end
fprintf(fid, 'id,name,freq_mode,time_mode,nfeatures,feature_idx\n');
for k = 1:nCfg
    fi = cfgs(k).feature_idx(:).';
    fprintf(fid, '%d,%s,%s,%s,%d,"%s"\n', ...
        k, cfgs(k).name, cfgs(k).freq_mode, cfgs(k).time_mode, ...
        numel(fi), strrep(mat2str(fi), '"', ''));
end
fclose(fid);

fprintf('Wrote %s (%d lines)\\n', out_txt, nCfg);
fprintf('Wrote %s\\n', out_csv);

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

% Fixed frequency setting: freqA only
freqA = [10 11];
freqName = 'freqA';

% Temporal feature settings to compare
timeSets = { ...
    [], ...          % no temporal feature
    17:20, ...       % prepOH only
    16, ...          % PT only
    [16 17:20] ...   % PT + prepOH
};

timeNames = { ...
    'none', ...
    'prepOH', ...
    'PT', ...
    'PT+prepOH'
};

cfgs = struct('name', {}, 'feature_idx', {}, 'freq_mode', {}, 'time_mode', {});

ic = 0;
for it = 1:numel(timeSets)
    ic = ic + 1;

    feat = [mandatory freqA timeSets{it}];
    feat = unique(feat, 'stable');

    cfgs(ic).freq_mode = freqName;
    cfgs(ic).time_mode = timeNames{it};
    cfgs(ic).feature_idx = feat;
    cfgs(ic).name = sprintf('%s__%s', freqName, timeNames{it});
end

end