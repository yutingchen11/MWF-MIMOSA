function cfg_id = resolve_cfg_id_loss(cfg_name, cfgs)
% English comments only

cfg_id = [];

if isnumeric(cfg_name) && isscalar(cfg_name)
    idx = double(cfg_name);
    if idx >= 1 && idx <= numel(cfgs) && abs(idx - round(idx)) < eps
        cfg_id = round(idx);
        return;
    end
end

cfg_name_str = char(string(cfg_name));
cfg_names = {cfgs.name};
cfg_id = find(strcmp(cfg_names, cfg_name_str), 1, 'first');

if isempty(cfg_id)
    idx = str2double(cfg_name_str);
    if ~isnan(idx) && isfinite(idx) && idx >= 1 && idx <= numel(cfgs) && abs(idx - round(idx)) < eps
        cfg_id = round(idx);
    end
end

if isempty(cfg_id)
    error('Unknown cfg_name: %s', cfg_name_str);
end
end