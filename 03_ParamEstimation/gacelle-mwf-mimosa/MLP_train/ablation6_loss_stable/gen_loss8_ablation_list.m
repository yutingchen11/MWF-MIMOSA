function gen_loss8_ablation_list(out_txt, out_csv)
% English comments only

if nargin < 1 || isempty(out_txt)
    out_txt = 'ablation_list.txt';
end
if nargin < 2 || isempty(out_csv)
    out_csv = 'ablation_list.csv';
end

cfgs = build_loss8_ablation_configs_local();
nCfg = numel(cfgs);

fid = fopen(out_txt, 'w');
if fid < 0
    error('Cannot open %s for writing.', out_txt);
end
for k = 1:nCfg
    fprintf(fid, '%s\n', cfgs(k).name);
end
fclose(fid);

fid = fopen(out_csv, 'w');
if fid < 0
    error('Cannot open %s for writing.', out_csv);
end

fprintf(fid, 'id,name,loss_abbr,residual_name,residual_mode,use_weight,use_norm,epsA,tau,nfeatures,feature_idx\n');
for k = 1:nCfg
    fi = cfgs(k).feature_idx(:).';
    fprintf(fid, '%d,%s,%s,%s,%d,%d,%d,%.8g,%.8g,%d,"%s"\n', ...
        k, cfgs(k).name, cfgs(k).loss_abbr, cfgs(k).residual_name, ...
        cfgs(k).residual_mode, cfgs(k).use_weight, cfgs(k).use_norm, ...
        cfgs(k).epsA, cfgs(k).tau, numel(fi), mat2str(fi));
end
fclose(fid);

fprintf('Wrote %s (%d lines)\n', out_txt, nCfg);
fprintf('Wrote %s\n', out_csv);

end

function cfgs = build_loss8_ablation_configs_local()
% English comments only

fixedFeatureIdx = [1:11 16:20];
tau  = 1e-2;
epsA = 1e-4;

defs = { ...
    'LOSS01__WRCL', 'WRCL', 'Charb', 1, 1, 1; ...
    'LOSS02__WCL',  'WCL',  'Charb', 1, 1, 0; ...
    'LOSS03__RCL',  'RCL',  'Charb', 1, 0, 1; ...
    'LOSS04__CL',   'CL',   'Charb', 1, 0, 0; ...
    'LOSS05__WRL1', 'WRL1', 'L1',    2, 1, 1; ...
    'LOSS06__WL1',  'WL1',  'L1',    2, 1, 0; ...
    'LOSS07__RL1',  'RL1',  'L1',    2, 0, 1; ...
    'LOSS08__L1',   'L1',   'L1',    2, 0, 0; ...
    };

cfgs = struct('name', {}, 'loss_abbr', {}, 'residual_name', {}, ...
              'residual_mode', {}, 'use_weight', {}, 'use_norm', {}, ...
              'tau', {}, 'epsA', {}, 'feature_idx', {});

for k = 1:size(defs, 1)
    cfgs(k).name          = defs{k,1}; %#ok<AGROW>
    cfgs(k).loss_abbr     = defs{k,2}; %#ok<AGROW>
    cfgs(k).residual_name = defs{k,3}; %#ok<AGROW>
    cfgs(k).residual_mode = defs{k,4}; %#ok<AGROW>
    cfgs(k).use_weight    = defs{k,5}; %#ok<AGROW>
    cfgs(k).use_norm      = defs{k,6}; %#ok<AGROW>
    cfgs(k).tau           = tau; %#ok<AGROW>
    cfgs(k).epsA          = epsA; %#ok<AGROW>
    cfgs(k).feature_idx   = fixedFeatureIdx; %#ok<AGROW>
end
end