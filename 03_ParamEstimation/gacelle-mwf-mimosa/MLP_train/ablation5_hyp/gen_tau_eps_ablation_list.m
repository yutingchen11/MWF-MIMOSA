function gen_tau_eps_ablation_list(out_txt, out_csv)
% English comments only

if nargin < 1 || isempty(out_txt)
    out_txt = 'ablation_list.txt';
end
if nargin < 2 || isempty(out_csv)
    out_csv = 'ablation_list.csv';
end

cfgs = build_tau_eps_ablation_configs_local();
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

fprintf(fid, 'id,name,loss_name,loss_mode,use_weight,use_norm,epsA,tau,nfeatures,feature_idx\n');
for k = 1:nCfg
    fi = cfgs(k).feature_idx(:).';
    fprintf(fid, '%d,%s,%s,%d,%d,%d,%.8g,%.8g,%d,"%s"\n', ...
        k, cfgs(k).name, cfgs(k).loss_name, cfgs(k).loss_mode, ...
        cfgs(k).use_weight, cfgs(k).use_norm, cfgs(k).epsA, cfgs(k).tau, ...
        numel(fi), mat2str(fi));
end
fclose(fid);

fprintf('Wrote %s (%d lines)\n', out_txt, nCfg);
fprintf('Wrote %s\n', out_csv);

end

function cfgs = build_tau_eps_ablation_configs_local()
% English comments only

fixedFeatureIdx = [1:11 16:20];

tauList = [1e-3, 3e-3, 1e-2, 3e-2, 1e-1];
epsList = [1e-6, 1e-5, 1e-4, 1e-3];

cfgs = struct('name', {}, 'feature_idx', {}, ...
              'loss_name', {}, 'loss_mode', {}, ...
              'use_weight', {}, 'use_norm', {}, ...
              'epsA', {}, 'tau', {});

ic = 0;
for itau = 1:numel(tauList)
    for ieps = 1:numel(epsList)
        ic = ic + 1;

        tauVal = tauList(itau);
        epsVal = epsList(ieps);

        cfgs(ic).feature_idx = fixedFeatureIdx; %#ok<AGROW>
        cfgs(ic).loss_name   = 'relCplx_full'; %#ok<AGROW>
        cfgs(ic).loss_mode   = 1; %#ok<AGROW>
        cfgs(ic).use_weight  = 1; %#ok<AGROW>
        cfgs(ic).use_norm    = 1; %#ok<AGROW>
        cfgs(ic).epsA        = epsVal; %#ok<AGROW>
        cfgs(ic).tau         = tauVal; %#ok<AGROW>

        cfgs(ic).name = sprintf('TAUEPS%02d__tau_%s__eps_%s', ...
            ic, val_to_tag_local(tauVal), val_to_tag_local(epsVal)); %#ok<AGROW>
    end
end
end

function tag = val_to_tag_local(x)
% English comments only

if x == 0
    tag = '0';
    return;
end

tag = sprintf('%.0e', x);
tag = strrep(tag, '+', '');
tag = strrep(tag, '-', 'm');
end