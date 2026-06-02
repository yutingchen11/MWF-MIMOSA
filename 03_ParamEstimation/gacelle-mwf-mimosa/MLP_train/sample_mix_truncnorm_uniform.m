function x = sample_mix_truncnorm_uniform(n, lo, hi, mu, sigma, p_gauss)
% English comments only
x = zeros(n,1,'single');

maskG = rand(n,1,'single') < single(p_gauss);
idxG  = find(maskG);
idxU  = find(~maskG);

% Uniform part
if ~isempty(idxU)
    x(idxU) = single(lo) + single(hi-lo) .* rand(numel(idxU),1,'single');
end

% Truncated normal via rejection sampling
if ~isempty(idxG)
    need = numel(idxG);
    out  = zeros(need,1,'single');
    got  = 0;
    while got < need
        m = max(ceil(2*(need-got)), 256);
        z = single(mu) + single(sigma) .* randn(m,1,'single');
        z = z(z >= lo & z <= hi);
        take = min(numel(z), need-got);
        if take > 0
            out(got+1:got+take) = z(1:take);
            got = got + take;
        end
    end
    x(idxG) = out;
end
end
