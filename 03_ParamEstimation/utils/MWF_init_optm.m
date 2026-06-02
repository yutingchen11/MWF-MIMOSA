%% ------------------------------------ fast approximate initialization for MWF

%--------------------------------------------------------------------------
% Basic settings
%--------------------------------------------------------------------------
num_acqs = param.ncontrast;

t1_entries = [150, 1000, 4500];   % MW, IEW, FW, in ms
t2_entries = [15, 70, 500];       % MW, IEW, FW, in ms

N = [size(img_all,1), size(img_all,2), size(img_all,3)];

%--------------------------------------------------------------------------
% Estimate R2' map
%--------------------------------------------------------------------------
R2p = 1 ./ T2s_map - 1 ./ T2_map;
R2p(R2p < 0) = 0;
R2p(~isfinite(R2p)) = 0;

%--------------------------------------------------------------------------
% Brain mask
%--------------------------------------------------------------------------
iMag = sqrt(sum(abs(img_all).^2, 4));

matrix_size = [size(iMag,1), size(iMag,2), size(iMag,3)];
voxel_size = [1 1 1];

iMag1 = zeros(matrix_size, 'like', iMag);
iMag1(60:end,:,:) = iMag(60:end,:,:);

mask_brain = BET(iMag1, matrix_size, voxel_size);
msk = mask_brain ~= 0;

save('mask_brain.mat', 'mask_brain');

%--------------------------------------------------------------------------
% Coarse R2' binning instead of using all unique T2* combinations
%--------------------------------------------------------------------------
%= R2p is in 1/ms if T2 maps are in ms. Convert to Hz.
R2p_Hz = R2p * 1e3;
R2p_Hz(~isfinite(R2p_Hz)) = 0;
R2p_Hz(R2p_Hz < 0) = 0;

%  =Median filtering stabilizes the initialization.
R2p_Hz_filt = medfilt3(R2p_Hz, [5 5 5]);

%  Coarse bins are enough for initialization.
r2p_bins_Hz = 0:2:30;

R2p_vec = R2p_Hz_filt(:);
R2p_vec(R2p_vec < min(r2p_bins_Hz)) = min(r2p_bins_Hz);
R2p_vec(R2p_vec > max(r2p_bins_Hz)) = max(r2p_bins_Hz);

dist_r2p = abs(bsxfun(@minus, R2p_vec, r2p_bins_Hz));
[~, t2s_bin_vec] = min(dist_r2p, [], 2);

t2s_bin_map = reshape(t2s_bin_vec, N);
t2s_bin_map(~msk) = 0;

%--------------------------------------------------------------------------
% Coarse B1 binning
%--------------------------------------------------------------------------
thre_high = 1.35;
thre_low  = 0.65;

img_b1_fast = img_b1;
img_b1_fast(~isfinite(img_b1_fast)) = 1;
img_b1_fast(img_b1_fast > thre_high) = thre_high;
img_b1_fast(img_b1_fast < thre_low)  = thre_low;
img_b1_fast(~msk) = 0;

%  : Coarse B1 grid for initialization.
% Use 0.65:0.05:1.35 if you want a closer result to the original code.
b1_val = 0.65:0.10:1.35;
num_b1 = length(b1_val);

b1_vec = img_b1_fast(:);
dist_b1 = abs(bsxfun(@minus, b1_vec, b1_val));
[~, b1_idx_vec] = min(dist_b1, [], 2);

b1_idx_map = reshape(b1_idx_vec, N);
b1_idx_map(~msk) = 0;

%--------------------------------------------------------------------------
% Build a coarse PV dictionary only once
%--------------------------------------------------------------------------
% Directly sample the simplex PV1 + PV2 + PV3 = 1.
% PV1: myelin water fraction
% PV2: intra/extracellular water fraction
% PV3: free water fraction

PV1_grid = 0.00:0.04:0.32;
PV3_grid = 0.00:0.05:0.70;

PV_dict = zeros(3, length(PV1_grid) * length(PV3_grid), 'single');

cnt = 0;
for pp = 1:length(PV1_grid)
    for rr = 1:length(PV3_grid)

        pv1 = PV1_grid(pp);
        pv3 = PV3_grid(rr);
        pv2 = 1 - pv1 - pv3;

        if pv2 >= 0.02 && pv2 <= 1
            cnt = cnt + 1;
            PV_dict(:, cnt) = single([pv1; pv2; pv3]);
        end
    end
end

PV_dict = PV_dict(:, 1:cnt);
pv = PV_dict.';
length_dict = size(pv, 1);

fprintf('Fast PV dictionary size: %d\n', length_dict);

%--------------------------------------------------------------------------
% Prepare image matrix: [num_acqs, Nvox]
%--------------------------------------------------------------------------
img_contrasts = abs(img_all);

if size(img_contrasts, 4) ~= num_acqs
    warning('size(img_all,4) is different from param.ncontrast.');
    num_acqs = size(img_contrasts, 4);
end

img_mat = reshape(img_contrasts, [], num_acqs);
img_mat = img_mat.';

%  Normalize each voxel signal for dictionary matching.
img_power = sum(abs(img_mat).^2, 1);
img_norm = sqrt(img_power);

valid_vox = img_norm > 0;
img_mat(:, valid_vox) = bsxfun(@rdivide, img_mat(:, valid_vox), img_norm(valid_vox));

%--------------------------------------------------------------------------
% Output maps
%--------------------------------------------------------------------------
PV1_map = zeros(N, 'single');
PV2_map = zeros(N, 'single');
PV3_map = zeros(N, 'single');

%  Chunk size avoids very large matrix multiplication.
chunk_size = 50000;

%--------------------------------------------------------------------------
% Dictionary matching over coarse T2* and B1 bins
%--------------------------------------------------------------------------
iniTime1 = clock;

for ind_t2sn = 1:length(r2p_bins_Hz)

    msk_t2s = msk & (t2s_bin_map == ind_t2sn);
    num_vox_t2s = sum(msk_t2s(:));

    if num_vox_t2s == 0
        continue;
    end

    %  : Convert the R2' bin back to 1/ms.
    r2p_bin_per_ms = r2p_bins_Hz(ind_t2sn) / 1e3;

    %  : Approximate component-specific T2* values.
    t2s_entries = round(1 ./ (r2p_bin_per_ms + 1 ./ t2_entries));
    t2s_entries(t2s_entries < 1) = 1;

    T1_entries = t1_entries(:);
    T2_entries = t2_entries(:);
    T2s_entries = t2s_entries(:);

    t1t2t2s_lut = cat(2, T1_entries, T2_entries, T2s_entries);

    %  : Remove physically invalid entries.
    keep_idx = t1t2t2s_lut(:,1) >= t1t2t2s_lut(:,2);
    keep_idx = keep_idx & (t1t2t2s_lut(:,2) >= t1t2t2s_lut(:,3));

    t1t2_lut_prune = t1t2t2s_lut(keep_idx, :);

    if isempty(t1t2_lut_prune)
        continue;
    end

    % fprintf('R2p bin %.1f Hz, component dictionary size: %d, voxels: %d\n', ...
    %     r2p_bins_Hz(ind_t2sn), size(t1t2_lut_prune,1), num_vox_t2s);

    %----------------------------------------------------------------------
    % Simulate component signals for this T2* bin
    %----------------------------------------------------------------------
    signal_conv_fit = zeros(size(t1t2_lut_prune,1), num_acqs, num_b1, 'single');

    inv_eff = [0.5 0.85 0.96];

    for bb = 1:num_b1

        [Mz_, Mxy_] = sim_mwf_MIMOSA_bp( ...
            TR, ...
            alpha_deg, ...
            esp, ...
            turbo_factor, ...
            t1t2_lut_prune(:,1) * 1e-3, ...
            t1t2_lut_prune(:,2) * 1e-3, ...
            num_reps, ...
            echo2use, ...
            TR_mte, ...
            esp_mte, ...
            TEs, ...
            t1t2_lut_prune(:,3) * 1e-3, ...
            gap_between_readouts, ...
            time2relax_at_the_end, ...
            b1_val(bb), ...
            inv_eff.' );

        temp_conv = abs(Mxy_(:,:,end).');

        for nn = 1:size(temp_conv,1)
            temp_power = sum(abs(temp_conv(nn,:)).^2);
            temp_norm = sqrt(temp_power);

            if temp_norm > 0
                temp_conv(nn,:) = temp_conv(nn,:) / temp_norm;
            end
        end

        signal_conv_fit(:,:,bb) = single(temp_conv);
    end

    %----------------------------------------------------------------------
    % Convert component dictionary to PV dictionary
    %----------------------------------------------------------------------
    dict_cnv = permute(signal_conv_fit, [2, 1, 3]);   % [num_acqs, Ncomp, Nb1]

    dict_pv_full = zeros(num_acqs, length_dict, num_b1, 'single');

    for bb = 1:num_b1

        dict_pv_bb = single(dict_cnv(:,:,bb)) * PV_dict;

        dict_power = sum(abs(dict_pv_bb).^2, 1);
        dict_norm = sqrt(dict_power);
        dict_norm(dict_norm == 0) = 1;

        dict_pv_bb = bsxfun(@rdivide, dict_pv_bb, dict_norm);

        dict_pv_full(:,:,bb) = dict_pv_bb;
    end

    dict_use = permute(dict_pv_full, [2, 1, 3]);   % [Ndict, num_acqs, Nb1]

    %----------------------------------------------------------------------
    % Match voxels in this T2* bin and each B1 bin
    %----------------------------------------------------------------------
    for bb = 1:num_b1

        msk_use = msk_t2s & (b1_idx_map == bb);
        vox_idx = find(msk_use);

        if isempty(vox_idx)
            continue;
        end

        num_vox = length(vox_idx);

        for st = 1:chunk_size:num_vox

            ed = min(st + chunk_size - 1, num_vox);
            vox_idx_chunk = vox_idx(st:ed);

            img_chunk = img_mat(:, vox_idx_chunk);

            res = dict_use(:,:,bb) * img_chunk;

            [~, max_idx] = max(abs(res), [], 1);

            PV1_map(vox_idx_chunk) = single(pv(max_idx, 1));
            PV2_map(vox_idx_chunk) = single(pv(max_idx, 2));
            PV3_map(vox_idx_chunk) = single(pv(max_idx, 3));
        end
    end
end

%--------------------------------------------------------------------------
% Optional visualization
%--------------------------------------------------------------------------
imagesc3d2(PV1_map, s(PV1_map)/2 + [35 5 -10], 201, [180,180,180], [0,0.3]);
setGcf(.5);
colormap gray;

imagesc3d2(PV2_map, s(PV2_map)/2 + [35 5 -10], 202, [180,180,180], [0,1]);
setGcf(.5);
colormap gray;

imagesc3d2(PV3_map, s(PV3_map)/2 + [35 5 -10], 203, [180,180,180], [0,1]);
setGcf(.5);
colormap gray;
% 
% save('MWF_initialization_fast.mat', 'PV1_map', 'PV2_map', 'PV3_map', ...
%     'r2p_bins_Hz', 'b1_val', 'PV_dict');