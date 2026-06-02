function gen_ablation_list(out_txt, out_csv)
% English comments only

if nargin < 1 || isempty(out_txt)
    out_txt = 'ablation_list.txt';
end
if nargin < 2 || isempty(out_csv)
    out_csv = 'ablation_list.csv';
end

cfgs = build_ablation_configs_arch_feat();
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

fprintf(fid, 'id,name,freq_mode,time_mode,neuron_div,numNeurons,nfeatures,feature_idx\n');
for k = 1:nCfg
    fi = cfgs(k).feature_idx(:).';
    nn = cfgs(k).numNeurons(:).';
    fprintf(fid, '%d,%s,%s,%s,%d,"%s",%d,"%s"\n', ...
        cfgs(k).id, cfgs(k).name, cfgs(k).freq_mode, cfgs(k).time_mode, ...
        cfgs(k).neuron_div, mat2str(nn), cfgs(k).nfeatures, mat2str(fi));
end
fclose(fid);

fprintf('Wrote %s (%d lines)\n', out_txt, nCfg);
fprintf('Wrote %s\n', out_csv);

end