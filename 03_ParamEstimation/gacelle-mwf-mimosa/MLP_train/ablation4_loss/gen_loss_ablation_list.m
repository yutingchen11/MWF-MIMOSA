function gen_loss_ablation_list(out_txt)
% English comments only

if nargin < 1 || isempty(out_txt)
    out_txt = 'ablation_list.txt';
end

cfgs = build_ablation_configs_loss();

fid = fopen(out_txt, 'w');
if fid < 0
    error('Cannot open %s for writing.', out_txt);
end

for k = 1:numel(cfgs)
    fprintf(fid, '%s\n', cfgs(k).name);
end

fclose(fid);

fprintf('Wrote %s (%d lines)\n', out_txt, numel(cfgs));
end