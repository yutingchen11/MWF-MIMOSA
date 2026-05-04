% generate dict first by run gen_MWF_MIMOSA_dict_1mm.m
% install sepia and 
% χ-separation Tool: https://github.com/SNU-LIST/chi-separation
 addpath(genpath('/gacelle-mwf-mimosa'));   % https://github.com/mriphysics/EPG-X

%% load dict
dict = h5read('dict/dict_mwf_mimosa_bp_50b1_126T2s_double.h5','/dict');
t2s_entries  = [1:1:100, 105:5:200 210:50:500];
%% load look up table
addpath(genpath('utils'))
load('dict/ielookup_4mimosa.mat');
%% load recon data
load('../Recon/zsssl_recon_1mm/output/mimosa_R8_zsssl_tfv2_epoch20000.mat')
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% set protocol  
%--------------------------------------------------------------------------

param.esp             = 5.8 * 1e-3;
param.turbo_factor    = 127;

param.TR          = 4500e-3 - 5.8*127e-3 + 27.7*127e-3 - 162.4*2e-3;
param.alpha_deg   = 4;
param.num_reps    = 5;
param.echo2use    = 1;
param.gap_between_readouts    = 900e-3;
param.time2relax_at_the_end   = 0;



% ##### mgre
param.TE_mte = [1.98:2.58:26.5].*1e-3;
param.TR_mte = 27.7e-3;
param.esp_mte = 2.58e-3;
nechoes = length(param.TE_mte);
param.ncontrast =4+nechoes;%nacq

TR                      = param.TR;
num_reps                = 5;
echo2use                = 1;
gap_between_readouts    = 900e-3;
time2relax_at_the_end   = 0;
alpha_deg = 4;
esp             = param.esp;
turbo_factor    = param.turbo_factor;

TR_mte = param.TR_mte;
esp_mte = param.esp_mte;
TEs = param.TE_mte;
nechoes = length(TEs);

img_all = zeros([240,224,192,14]);
img_all(50:230,:,:,:) = img_zsssl;% use when skip neck; use 1:240 if need neck
T1_all = zeros(size(img_all,1),size(img_all,2),size(img_all,3));
T2_all = zeros(size(img_all,1),size(img_all,2),size(img_all,3));
PD_all = zeros(size(img_all,1),size(img_all,2),size(img_all,3));
T2s_all = zeros(size(img_all,1),size(img_all,2),size(img_all,3));
IE_all = zeros(size(img_all,1),size(img_all,2),size(img_all,3));


img = abs(img_all);
%--------------------------------------------------------------------------
%% Load b1 map and resize to qalas, for spetial reduced B1 slices, 42, for zju prisma B1map product
%--------------------------------------------------------------------------

matrix_size = [size(img,1),size(img,2),size(img,3)];

img_b1_load = niftiread('dcm_b1_tfl_b1map_20260302130524_6001.nii');
img_b1 = double(img_b1_load)/800;       % reference fa 800
% img_b1 = flip(img_b1, 1);  
img_b1 = permute(img_b1, [2 1 3]); 
img_b1 = flip(img_b1,2); %  

img_b1 = imresize3(img_b1, [240,240,195]); %  to match img size
img_b1 = img_b1(:,9:end-8,2:end-2);

imagesc3d2(img_b1, s(img_b1)/2, 51, [-0,-0,0], [0,2])
imagesc3d2(img(:,:,:,1), s(img(:,:,:,1))/2, 52, [-0,-0,0], [0,2e-4])

save('B1_4mimosa.mat','img_b1')
%% Select slice
%--------------------------------------------------------------------------

N = [size(img, 1) size(img, 2) size(img, 3)];

% % ROI mask for Dicoms
for slc_select = 1:N(3)
   msk(:,:,slc_select) = imfill(rsos(img(:,:,slc_select,:),4) > 0, 'holes'); % or 50(invivo) / 100(nist) thre
end

%--------------------------------------------------------------------------
%% threshold high and low B1 values: use raw b1 map without polyfit
%--------------------------------------------------------------------------

thre_high   = 1.35;
thre_low    = 0.65;

temp        = img_b1 .* msk;

temp(temp > thre_high)  = thre_high;
temp(temp < thre_low)   = thre_low;

temp        = temp .* msk;
img_b1      = temp .* msk;


%--------------------------------------------------------------------------
%% create masks for each b1 value
%--------------------------------------------------------------------------

num_b1_bins = 50; % default: 50

b1_val      = linspace( min(img_b1(msk==1)), max(img_b1(msk==1)), num_b1_bins );

sum_msk     = sum(msk(:));

if length(b1_val) == 1
    % do not use b1 correction
    msk_b1 = msk;
else
    
    msk_b1 = zeross([N,length(b1_val)]);
    
    for t = 1:length(b1_val)
        if t > 1
            msk_b1(:,:,:,t) = (img_b1 <= b1_val(t)) .* (img_b1 > b1_val(t-1));
        else
            msk_b1(:,:,:,t) = msk.*(img_b1 <= b1_val(t));
        end
        
        percent_bin = sum(sum(sum(msk_b1(:,:,:,t),1),2),3) / sum_msk;
        
        if t == length(b1_val)
            msk_b1(:,:,:,t) = img_b1 > b1_val(t-1);
        end
        
        msk_b1(:,:,:,t) = msk_b1(:,:,:,t) .* msk;
    end
end


%% estimate T2* separatly using vapro

image_series = abs(squeeze(img(:,:,:,5:end)));
TE = TEs*1e3;
[Np, Nf, Nfr, Nte] = size(image_series);
Nsize = [Np, Nf, Nfr];
possibleT2Values = [1:500]; % the range of possible T2 values
[T2std, Mstd] = t2EstiVAPRO(reshape(image_series, Np*Nf*Nfr, [])', TE, possibleT2Values); 
Mstd  = reshape(Mstd, Np, Nf,Nfr);
T2std = reshape(T2std, Np, Nf, Nfr);

imagesc3d2(sq(T2std), s(sq(T2std))/2+[30 0 0], 99, [180 180 180], [0 200]);colormap hot

%% distribution of T2* 

msk = T2std>1;
A_flat = T2std(msk);
data = A_flat(:);
values = unique(data); %
counts = histc(data(:), values); % 

figure;
plot(values, counts);

%% 
%--------------------------------------------------------------------------
% threshold high and low T2* values: use T2* map
%--------------------------------------------------------------------------

thre_high   = max(values);
thre_low    = min(values);

temp        = T2std .* msk;

temp(temp > thre_high)  = thre_high;
temp(temp < thre_low)   = thre_low;

temp        = temp .* msk;
T2std      = temp .* msk;


%--------------------------------------------------------------------------
% create masks for each T2s value
%--------------------------------------------------------------------------


T2s_val      = values;
%b1_val      = 1; % no b1 correction
sum_msk     = sum(msk(:));

if length(T2s_val) == 1
    % do not use b1 correction
    msk_t2s = msk;
else
    
    msk_t2s = zeross([N,length(T2s_val)]);
    
    for t = 1:length(T2s_val)
        if t > 1
            msk_t2s(:,:,:,t) = (T2std <= T2s_val(t)) .* (T2std > T2s_val(t-1));
        else
            msk_t2s(:,:,:,t) = msk.*(T2std <= T2s_val(t));
        end
        
        percent_bin = sum(sum(sum(msk_t2s(:,:,:,t),1),2),3) / sum_msk;
        
        if t == length(T2s_val)
            msk_t2s(:,:,:,t) = T2std > T2s_val(t-1);
        end
        
        msk_t2s(:,:,:,t) = msk_t2s(:,:,:,t) .* msk;
    end
end

%--------------------------------------------------------------------------
%% create look up table 
%--------------------------------------------------------------------------
tic
clc
fprintf('generating dictionaries...\n');

t1_entries  = [5:10:3000, 3100:50:5000];
t2_entries  = [1:2:350, 370:20:1000, 1100:100:3000];


T1_entries  = repmat(t1_entries.', [1,length(t2_entries)]).';
T1_entries  = T1_entries(:);
  
T2_entries  = repmat(t2_entries.', [1,length(t1_entries)]);
T2_entries  = T2_entries(:);


t1t2_lut    = cat(2, T1_entries, T2_entries);

% remove cases where T2>T1
idx = 0;
for t = 1:length(t1t2_lut)
    if t1t2_lut(t,1) < t1t2_lut(t,2)
        idx = idx+1;
    end
end

t1t2_lut_prune = zeross( [length(t1t2_lut) - idx, 2] );


ie_lkp = reshape((ielookup.ies_mtx),[length(t1t2_lut),1]);
ie_lkp_prune = zeross( [length(t1t2_lut) - idx, 1] );
%% prune
idx = 0;
for t = 1:length(t1t2_lut)
    if t1t2_lut(t,1) >= t1t2_lut(t,2)
        idx = idx+1;
        t1t2_lut_prune(idx,:) = t1t2_lut(t,:);
        ie_lkp_prune(idx,:) = ie_lkp(t,:);
    end
end

disp(['dictionary entries: ', num2str(length(t1t2_lut_prune))])

%%
%--------------------------------------------------------------------------
% dictionary fit -> in each slice, bin voxels based on b1 value
%--------------------------------------------------------------------------
% img = img_sim;
length_dict = length(t1t2_lut_prune);
fprintf('fitting...\n');

estimate_pd_map = 1;     % set to 1 to estiamte PD map (makes it slower)

T1_map = zeross(N);
T2_map = zeross(N);
% T2s_map = zeross(N);

PD_map = zeross(N);     % proton density-> won't be estimated if estimate_pd_map is set to 0
IE_map = zeross(N);     % inversion efficiency

iniTime         = clock;
iniTime0        = iniTime;


% parallel computing
% delete(gcp('nocreate'))
% c = parcluster('local');
% total_cores = c.NumWorkers;
% parpool(min(ceil(total_cores*0.5), 15))

for slc_select = 1:N(3)
    disp(num2str(slc_select))

    tic
    % for b = 1:length(b1_val)
    for t2ss = 1:length(T2s_val)
        msk_slc_t2s = msk_t2s(:,:,slc_select,t2ss);
        % msk_slc = new_mask(:,:,b);
        num_vox = sum(msk_slc_t2s(:)~=0);

        if num_vox > 0

            for b = 1:length(b1_val)
                msk_slc_b1 = msk_b1(:,:,slc_select,b);
                num_vox_b1 = sum(msk_slc_b1(:)~=0);
                msk_slc = msk_slc_t2s.*msk_slc_b1;
                % te = msk_slc + te;
                num_vox_all = sum(msk_slc(:)~=0);
                if num_vox_all > 0
                    
                    img_slc = zeros([14, num_vox_all]);
    
                    for t = 1:14
                        temp = sq(img(:,:,slc_select,t));
                        img_slc(t,:) = temp(msk_slc~=0);
                    end
          


                        [~, t2s_idx] = min(abs(t2s_entries - T2s_val(t2ss)));
                        if b1_val == 1
                            res = dict(:,:,t2s_idx,25) * img_slc; % dot product
                        else
                            res = dict(:,:,t2s_idx,b) * img_slc;
                        end
                        
        
                        % find T1, T2 values
                        [~, max_idx] = max(abs(res), [], 1);
        
                        max_idx_t1t2 = mod(max_idx, length_dict);
                        max_idx_t1t2(max_idx_t1t2==0) = length_dict;
        
                        res_map = t1t2_lut_prune(max_idx_t1t2,:);
        
 
                        ie_to_use = ie_lkp_prune(max_idx_t1t2,1);
 
        
                        if estimate_pd_map
                            [Mz_sim, Mxy_sim] = sim_mwf_MIMOSA_bp(TR, alpha_deg, esp, turbo_factor, res_map(:,1)*1e-3, res_map(:,2)*1e-3, num_reps, echo2use,TR_mte,esp_mte,TEs ,T2s_val(t2ss)*1e-3, gap_between_readouts, time2relax_at_the_end, b1_val(b), ie_to_use);
                           
                        end
                        
        
                    t1_map = zeross(N(1:2));
                    t1_map(msk_slc==1) = res_map(:,1);
        
                    t2_map = zeross(N(1:2));
                    t2_map(msk_slc==1) = res_map(:,2);
                    
        
                    ie_map = zeross(N(1:2));
                    ie_map(msk_slc==1) = ie_to_use;
        
                    if estimate_pd_map
                        Mxy_sim_use = abs(Mxy_sim(:,:,end));
        
                        scl = zeross([num_vox_all,1]);
        
                        for idx = 1:size(Mxy_sim_use,2)
                            scl(idx) = Mxy_sim_use(:,idx) \ img_slc(:,idx);
                        end
        
                        pd_map = zeross(N(1:2));
                        pd_map(msk_slc~=0) = scl;
                        PD_map(:,:,slc_select) = PD_map(:,:,slc_select) + pd_map;
                    end
        
                    T1_map(:,:,slc_select) = T1_map(:,:,slc_select) + t1_map;
                    T2_map(:,:,slc_select) = T2_map(:,:,slc_select) + t2_map;
                    IE_map(:,:,slc_select) = IE_map(:,:,slc_select) + ie_map;
    
                end
            end
        end
    end
    toc
end
%% 3D 
T2s_map = T2std;
imagesc3d2(sq(T1_map), s(sq(T1_map(:,:,:)))/2+[30 0 0], 1, [180 180 180], [0 3000]);colormap hot
imagesc3d2(sq(T2_map), s(sq(T1_map(:,:,:)))/2+[30 0 0], 2, [180 180 180], [0 200]);colormap hot
imagesc3d2(sq(T2s_map), s(sq(T1_map(:,:,:)))/2+[30 0 0], 3, [180 180 180], [0 200]);colormap hot
imagesc3d2(sq(PD_map), s(sq(T1_map(:,:,:)))/2+[30 0 10], 4, [180 180 180], [0 max(abs(PD_map(:)))]);colormap gray
imagesc3d2(sq(IE_map), s(sq(T1_map(:,:,:)))/2+[30 0 10], 5, [180 180 180], [0 1]);colormap gray


imagesc3d2(sq(fliplr(permute(T1_map,[3 2 1]))), s(sq(permute(T1_map,[2 1 3])))/2+[0 0 48], 1, [90 90 -90], [0 3000]);colormap hot
imagesc3d2(sq(fliplr(permute(T2_map,[3 2 1]))), s(sq(permute(T2_map,[2 1 3])))/2+[0 0 48], 2, [90 90 -90], [0 200]);colormap hot

%%
save('mapping_mwf_mimosa_R8_bp_zsssl.mat','T1_map','T2_map','T2s_map','PD_map','IE_map','img_all','-v7.3')
%% save nii files
is_save_nii = 1;

if is_save_nii
    nii_T1 = make_nii(sq((T1_map)));
    save_nii(nii_T1,'mimosa_sub5_R8_T1.nii');
    nii_T2 = make_nii(sq((T2_map)));
    save_nii(nii_T2,'mimosa_sub5_R8_T2.nii');
    nii_T2s = make_nii(sq((T2s_map)));
    save_nii(nii_T2s,'mimosa_sub5_R8_T2s.nii');
    nii_pd = make_nii(sq((PD_map)));
    save_nii(nii_pd,'mimosa_sub5_R8_pd.nii');
end

%% ------------------------------------initialization for MWF
R2p = (1/T2s_map - 1/T2_map);
R2p(R2p<0) = 0; 
R2p(isnan(R2p)) = 0; 
imagesc3d2(R2p*1e3, s(R2p)/2+[28 5 -10], 100, [0,180,180], [0,30]), setGcf(.5);colormap hot
num_acqs = param.ncontrast;
%% generating dictionaries for conventional and subspace mapping
% v13
t1_entries  = [150,1000,4500]; %   MW WM CSF
t2_entries  = [15,70,500]; % MW WM CSF
% t2s_entries  = [10,70,500];% WM/GM/CSF
%% mask
iMag = sqrt(sum(abs(img_all).^2,4));%echo time diemens sum
matrix_size = [size(iMag,1),size(iMag,2),size(iMag,3)];
voxel_size = [1 1 1];
iMag1 = zeros(matrix_size);
iMag1(60:end,:,:) = iMag(60:end,:,:);

mask_brain = BET(iMag1,matrix_size,voxel_size);
imagesc3d2(mask_brain.*iMag, s(sq(mask_brain))/2+[35 5 -10], 8, [90,90,90], [0 2e-3]);colormap gray
save('mask_brain.mat','mask_brain');
%% estimate T2s of 3 components
load('mask_brain.mat')
msk = mask_brain;

T2sn = round(1/(repmat(R2p.*mask_brain,[1 1 1 3]) + reshape(1./(t2_entries),[1 1 1 3]))).*msk;
% to save time
for ii = 1:3
    T2sn(:,:,:,ii) = medfilt3(sq(T2sn(:,:,:,ii)),[5 5 5]);
end

imagesc3d2(T2sn(:,:,:,1), s(R2p)/2+[35 5 -10], 101, [180,180,180], [0,20]), setGcf(.5);colormap hot
imagesc3d2(T2sn(:,:,:,2), s(R2p)/2+[35 5 -10], 102, [180,180,180], [0,100]), setGcf(.5);colormap hot
imagesc3d2(T2sn(:,:,:,3), s(R2p)/2+[35 5 -10], 103, [180,180,180], [0,500]), setGcf(.5);colormap hot

all_combinations = reshape(T2sn, [], 3);  

[T2sn_comb, ~, ic] = unique(all_combinations, 'rows');
 counts = histc(ic, 1:size(T2sn_comb,1));
 figure;plot(counts)
 figure;plot(counts(2:end-1))
%% bin combinations
%--------------------------------------------------------------------------
% create masks for each T2s value
%--------------------------------------------------------------------------
[Nx,Ny,Nz] = size(T2s_map);
N  = [size(img_all,1),size(img_all,2),size(img_all,3)];

if length(T2sn_comb) == 1
    % do not use b1 correction
    msk_t2sn = msk;
else
    
    msk_t2sn = zeross([N,length(T2sn_comb)]);
    
    for t = 1:length(T2sn_comb)
        current_combo = T2sn_comb(t, :);
        is_this_combo = (ic == t);  
        msk_t2sn(:,:,:,t) = reshape(is_this_combo, Nx, Ny, Nz);

        msk_t2sn(:,:,:,t) = msk_t2sn(:,:,:,t).*msk;
    end
end


PV1_map = zeros([Nx,Ny,Nz]);
PV2_map = zeros([Nx,Ny,Nz]);
PV3_map = zeros([Nx,Ny,Nz]);
%%        
for ind_t2sn = 1:length(T2sn_comb)
num_vox = sum(msk_t2sn(:,:,:,ind_t2sn)~=0,'all');
    if num_vox>0
        t2s_entries = T2sn_comb(ind_t2sn,:);
        
        inv_eff     = [0.5 0.85 0.96]; % 0.5:0.05:1.0
        inv_eff_sub = [0.5 0.85 0.96]; % 0.5:0.05:1.0
        b1_val      =  [0.65:0.05:1.35]; % 
        b1_val_sub  =  [0.65:0.05:1.35]; % 
        
        T1_entries  = t1_entries; % PV
        T1_entries  = T1_entries(:);
        
        T2_entries  = t2_entries; % PV
        T2_entries  = T2_entries(:);
        
        T2s_entries  = t2s_entries; % PV
        T2s_entries  = T2s_entries(:);
        
        t1t2t2s_lut    = cat(2, T1_entries, T2_entries, T2s_entries);
        
        % remove cases where T2>T1
        idx = 0;
        for t = 1:length(t1t2t2s_lut)
            if t1t2t2s_lut(t,1) < t1t2t2s_lut(t,2) || t1t2t2s_lut(t,2) < t1t2t2s_lut(t,3)
                idx = idx+1;
            end
        end
        
        t1t2_lut_prune = zeross([length(t1t2t2s_lut) - idx, 3]);
        
        idx = 0;
        for t = 1:length(t1t2t2s_lut)
            if t1t2t2s_lut(t,1) >= t1t2t2s_lut(t,2)
                idx = idx+1;
                t1t2_lut_prune(idx,:) = t1t2t2s_lut(t,:);
            end 
        end
        
        fprintf('dictionary size: %d \n', length(t1t2_lut_prune)*length(b1_val));
        
        E = turbo_factor * num_acqs;
        
        
        signal_conv_fit     = zeross([length(t1t2_lut_prune), num_acqs, length(b1_val)]);

        length_b1_val       = length(b1_val);
        length_b1_val_sub   = length(b1_val_sub);
        length_inv_eff      = length(inv_eff);
        length_inv_eff_sub  = length(inv_eff_sub);
        
        
        % for fitting
        cnt         = 0;
        iniTime1    = clock;
        
        for b1 = 1:length_b1_val
            % for ie = 1:length_inv_eff
                cnt             = cnt + 1;
                [Mz_, Mxy_]     = sim_mwf_MIMOSA_bp(TR, alpha_deg, esp, turbo_factor, t1t2_lut_prune(:,1)*1e-3, t1t2_lut_prune(:,2)*1e-3, num_reps, echo2use,TR_mte,esp_mte,TEs,  t1t2_lut_prune(:,3)*1e-3, gap_between_readouts, time2relax_at_the_end, b1_val(b1), inv_eff.');
                
                temp_conv       = abs(Mxy_(:,:,end).');
                
                for n = 1:size(temp_conv,1)
                    temp_conv(n,:)      = temp_conv(n,:) / sum(abs(temp_conv(n,:)).^2)^0.5;
                end
                
                signal_conv_fit(:,:,b1)  = temp_conv;
        
            % end
        end
        
        delete(gcp('nocreate'))
        fprintf('total elapsed time: %.1f sec\n\n',etime(clock,iniTime1));
        

        % for conventional fitting
        length_dict = length(t1t2_lut_prune);
        dict_cnv    = zeross([length_dict * 1 , num_acqs, length(b1_val)]);
        % dict_cnv = reshape(permute(signal_conv_fit,[1,3,4,5,2]),[length_dict * length(inv_eff) * length_t2s_val * length(b1_val),num_acqs]);
        for t = 1:1
            dict_cnv(1 + (t-1)*length_dict : t*length_dict, :, :) = signal_conv_fit(:,:,:,t);
        end

        
        %% PV dictionary-based
        
        dict_cnv_    = permute(dict_cnv,[2,1,3]);

        PV_dict_1   = 0.02:0.02:0.3; % myelin water
        PV_dict_2   = 0.02:0.02:1.0; % intra/extracellular water
        PV_dict_3   = 0.05:0.05:1.0; % free water
        
        PV_dict     = zeros(3,length(PV_dict_1)*length(PV_dict_2)*length(PV_dict_3));
        
        cnt = 0;
        
        for pp = 1:length(PV_dict_1)
            for qq = 1:length(PV_dict_2)
                for rr = 1:length(PV_dict_3)
                    cnt = cnt + 1;
                        PV_dict(:,cnt) = [PV_dict_1(pp),PV_dict_2(qq),PV_dict_3(rr)]./ ...
                            (PV_dict_1(pp)+PV_dict_2(qq)+PV_dict_3(rr));
                end
            end
        end
        PV_dict = unique(PV_dict','rows')';
        dict_pv_full=zeros(num_acqs, size(PV_dict,2), length(b1_val));
        for bb=1:length(b1_val)
            dict_pv_bb = squeeze(dict_cnv_(:,:,bb)) * PV_dict;
            dict_pv_full(:,:,bb) = dict_pv_bb ;
        end
        
        for cc = 1:size(dict_pv_full,2)
            % dict_pv_full(:,cc) = dict_pv_full(:,cc) / sum(abs(dict_pv_full(:,cc)).^2)^0.5;
            for bb = 1:length(b1_val)
                dict_pv_full(:,cc,bb) = dict_pv_full(:,cc,bb) / sum(abs(dict_pv_full(:,cc,bb)).^2)^0.5;
            end
        end
        %% PV dictionary-based, using all imgs
        
        img_contrasts = abs(img_all); % Nx Ny Nz echo 
        %%
        %--------------------------------------------------------------------------
        %% threshold high and low B1 values: use raw b1 map without polyfit
        %--------------------------------------------------------------------------
        
        thre_high   = 1.35;
        thre_low    = 0.65;
        
        temp        = img_b1 .* msk;
        
        temp(temp > thre_high)  = thre_high;
        temp(temp < thre_low)   = thre_low;
        
        temp        = temp .* msk;
        img_b1      = temp .* msk;
        
        
        %--------------------------------------------------------------------------
        %% create masks for each b1 value
        %--------------------------------------------------------------------------
        
        num_b1_bins = 15; % default: 50
        
        b1_val      = linspace( min(img_b1(msk==1)), max(img_b1(msk==1)), num_b1_bins );
        %b1_val      = 1; % no b1 correction
        sum_msk     = sum(msk(:));
        
        if length(b1_val) == 1
            % do not use b1 correction
            msk_b1 = msk;
        else
            
            msk_b1 = zeross([N,length(b1_val)]);
            
            for t = 1:length(b1_val)
                if t > 1
                    msk_b1(:,:,:,t) = (img_b1 <= b1_val(t)) .* (img_b1 > b1_val(t-1));
                else
                    msk_b1(:,:,:,t) = msk.*(img_b1 <= b1_val(t));
                end
                
                percent_bin = sum(sum(sum(msk_b1(:,:,:,t),1),2),3) / sum_msk;
                
                if t == length(b1_val)
                    msk_b1(:,:,:,t) = img_b1 > b1_val(t-1);
                end
                
                msk_b1(:,:,:,t) = msk_b1(:,:,:,t) .* msk;
            end
        end

        N  = [size(img_contrasts,1),size(img_contrasts,2),size(img_contrasts,3)];

        dict_use = permute(dict_pv_full,[2,1,3]);
        pv = permute(PV_dict,[2,1]);
        length_dict = length(pv);

        for ss = 1:size(img_contrasts,3)
            msk_slc_t2s = msk_t2sn(:,:,ss,ind_t2sn);
            num_vox_all = sum(msk_slc_t2s(:)~=0);
            if num_vox_all > 0
              for b = 1:length(b1_val)
                msk_slc_b1 = msk_b1(:,:,ss,b);
                num_vox_b1 = sum(msk_slc_b1(:)~=0);
                msk_slc = msk_slc_t2s.*msk_slc_b1;
                % te = msk_slc + te;
                num_vox_all = sum(msk_slc(:)~=0);
                if num_vox_all > 0


                    img_slc = zeros([14,num_vox_all]);
                    for t = 1:14
                        temp = sq(img_contrasts(:,:,ss,t));
                        img_slc(t,:) =  temp(msk_slc~=0);
                    end                    

                    res = dict_use(:,:,b) * img_slc;
                    [~, max_idx]    = max(abs(res), [], 1);
                    
                    max_idx_pv      = mod(max_idx, length_dict);
                    max_idx_pv(max_idx_pv==0) = length_dict;
                    
                    res_map         = pv(max_idx_pv,:);
                    
                    pv1_map = zeross(N(1:2));
                    pv1_map(msk_slc==1) = res_map(:,1);
                    
                    pv2_map = zeross(N(1:2));
                    pv2_map(msk_slc==1) = res_map(:,2);
        
                    pv3_map = zeross(N(1:2));
                    pv3_map(msk_slc==1) = res_map(:,3);
    
                % pv4_map = zeross(N(1:2));
                % pv4_map(msk_slc==1) = res_map(:,4);
                
                PV1_map(:,:,ss) = PV1_map(:,:,ss) + pv1_map;
                PV2_map(:,:,ss) = PV2_map(:,:,ss) + pv2_map;
                PV3_map(:,:,ss) = PV3_map(:,:,ss) + pv3_map;
            end
            fprintf('done\n')
              end
            end
        % end
        end
    end
end
imagesc3d2( PV1_map(:,:,:), s(PV1_map)/2+[30 0 0], 1, [180,180,180], [0,0.3]), setGcf(.5)

imagesc3d2( PV2_map(:,:,:), s(PV1_map)/2+[30 0 0], 2, [180,180,180], [0,1]), setGcf(.5)

imagesc3d2( PV3_map(:,:,:), s(PV1_map)/2+[30 0 0], 3, [180,180,180], [0,1]), setGcf(.5)


figure(10);
tiledlayout(1,3, 'TileSpacing', 'compact');
nexttile;imshow(sq(PV1_map(:,:,end/2-30)),[0,0.5]);title('Short T2 component');
nexttile;imshow(sq(PV2_map(:,:,end/2-30)),[0,1]);title('Intermediate T2 component');
nexttile;imshow(sq(PV3_map(:,:,end/2-30)),[0,1]);title('Long T2 component');


save('PV_fitting_mwf_mimosa_R8.mat','PV1_map','PV2_map','PV3_map')
%% release mem
clear dict
%--------------------------------------------------------------------------
%% GACELLE based MWF mapping
%--------------------------------------------------------------------------

TR = param.TR_mte;            % repetition time, in second
TE = param.TE_mte.'; % all echo times, in second, order should be the same as the way image data ared loaded

b1 = img_b1;
%% correct of phase offset
mosaic( angle(sq(img_all(150,:,:,end-9:end))), 2, 5, 101, '',[-pi pi],0);colormap jet %

algorParam.unwrap.isBipolarCorr = true;             
algorParam.unwrap.bipolarCorr.readoutDir = [1 0 0]; 
algorParam.unwrap.bipolarCorr.method = 'linear';  

algorParam = [];
algorParam.unwrap.unwrapMethod = 'ROMEO';

[bipolarCorr,FIT3D] = BipolarEddyCorrect(img_all(:,:,:,end-9:end),logical((mask_brain)),algorParam);

mosaic( angle(sq(bipolarCorr(150,:,:,:))), 2, 5, 102, '',[-pi pi],0);colormap jet %

img_cor = cat(4,img_all(:,:,:,1:4),bipolarCorr);
%% field estimation
cplx_img = img_cor(:,:,:,end-9:end);
[iFreq_raw, N_std ,relres ,p0] = Fit_ppm_complex_TE(cplx_img,TE);
imagesc3d2(sq(iFreq_raw(:,:,:)/2/pi/(TE(2)-TE(1))), s(sq(iFreq_raw(:,:,:,1,1)))/2+[35 0 0], 7, [180,180,180], [-100 100]);colormap parula% mosaic(sq(p0(:,:,:)/(TE(2)-TE(1))), 1,1,11, '', [-pi pi],180);colormap parula

totalField = iFreq_raw/2/pi/(TE(2)-TE(1));
pini = angle(cplx_img(:,:,:,1) ./ exp(1i*2*pi*totalField(:,:,:)*TE(1)));
imagesc3d2(pini, s(sq(pini))/2+[30 0 0], 8, [180,180,180], [-pi pi]);colormap parula
%% fixed parameters
kappa_mw                = 0.36; % Jung, NI., myelin water density
kappa_iew               = 0.86; % Jung, NI., intra-/extra-axonal water density
fixed_params.B0     	= 2.89;    % field strength, in tesla
fixed_params.rho_mw    	= kappa_mw/kappa_iew; % relative myelin water density
fixed_params.E      	= 0.02; % exchange effect in signal phase, in ppm
fixed_params.x_i      	= -0.1; % myelin isotropic susceptibility, in ppm
fixed_params.x_a      	= -0.1; % myelin anisotropic susceptibility, in ppm
fixed_params.B0dir      = [0 0 1];
%%
single_slc_check = 1;

if single_slc_check
slc = 148;

img = img_cor(slc,:,:,:);
mask_slc = mask_brain(slc,:,:);

y = img;
mask = mask_slc;

% itinital values
totalField_slc = totalField(slc,:,:);
pini_slc = pini(slc,:,:);
b1_slc = b1(slc,:,:);
scaleFactor = max(PD_map(:));
% PD_slc = PD_map(slc,:,:)/scaleFactor;
PD_slc = PD_map(slc,:,:);
PV1_slc = PV1_map(slc,:,:);
PV2_slc = PV2_map(slc,:,:);

pars0 = [];
pars0.fieldmap   = totalField_slc;   % Hz    % it is not neccessary to provide the total field for initial starting point, yet it can provide the fitting results
pars0.pini       = pini_slc;% rad
pars0.b1map     = b1_slc;
% [pars0.PDmap, ~] =  xy2voxels(PD_slc,mask_slc)/scaleFactor;
pars0.fmw_init = PV1_slc;
pars0.fiew_init = PV2_slc;

param.TE_flash = 2.29e-3;
param.TEs = param.TE_mte;
param.nechoes = numel(param.TE_mte);
%%
objGPU                     = gpuMIMOSAMWI_A05(TE.',fixed_params);


fitting = [];
fitting.param =  param;
fitting.iteration           = 2000;
fitting.initialLearnRate    = 0.01;% 0.01 ; 0.003
fitting.decayRate           = 0;
fitting.convergenceValue    = 1e-6;% slope of tol loss
fitting.tol                 = 1e-6;% tot loss, noise level; 
fitting.isDisplay           = 1;
fitting.weightPower         = 0.5;
fitting.patience          = 50;     % allow longer plateau before stopping; sometimes fitting of all voxels are not finished yet
% fitting.convergenceWindow = 80;     % larger window when smaller convengen => less sensitive to tiny oscillation
% fitting.decayRate         = 1e-3;% to reduce the lr
% fitting.lambda = 1e-4;

fitting.isWeighted = false;
fitting.debug         = 1;
fitting.enableComplex = false;% 
fitting.start = 'prior_mimosa';% 
fitting.isComplex = 1;

fitting.DIMWI.isFitIWF      = 1;
fitting.DIMWI.isFitFreqMW   = 1;
fitting.DIMWI.isFitFreqIEW   = 1;
fitting.DIMWI.isFitR2sIEW    = 1;

fitting.lambda = 3e-4;
fitting.regmap = {'MWF'};
fitting.TVmode = '3D';

fitting.lossFunction = 'l1';   

extraData = [];
extraData.freqBKG   = pars0.fieldmap / (gpuMIMOSAMWI.gyro*fixed_params.B0); % in ppm
extraData.pini      =  pars0.pini;
extraData.IWF       = pars0.fiew_init;

scale = prctile(PD_slc(:),98);

extraData.PD       = PD_slc/scale;
extraData.MWF       = PV1_slc;
extraData.b1       = b1_slc;

fitting.usingANN = true;   % MLP
fitting.model= 'invivo';
[out_full]    = objGPU.estimate(y/scale, mask, extraData, fitting);
figure(77);imshow(imrotate(sq(out_full.final.freqMW*gpuMIMOSAMWI.gyro*fixed_params.B0),180),[-20 20]);colormap jet;title('Freq MW')
figure(33);imshow(imrotate(sq(out_full.final.MWF),180),[0 0.3]);colormap gray;

end
%% 3D
Nx = size(img_all,2);
Ny = size(img_all,3);
Nz = size(img_all,1);

results_3D = struct();
results_3D.S0 = zeros(Nz, Nx, Ny);
results_3D.MWF = zeros(Nz, Nx, Ny);
results_3D.IEW = zeros(Nz, Nx, Ny);
results_3D.FW = zeros(Nz, Nx, Ny);  % 1 - MWF - IEW
results_3D.r2sMW = zeros(Nz, Nx, Ny);
results_3D.r2sIEW = zeros(Nz, Nx, Ny);
results_3D.freqMW = zeros(Nz, Nx, Ny);
results_3D.freqIEW = zeros(Nz, Nx, Ny);
results_3D.dfreqBKG = zeros(Nz, Nx, Ny);
results_3D.dpini = zeros(Nz, Nx, Ny);
results_3D.fieldmap_total = zeros(Nz, Nx, Ny);
results_3D.pini_total = zeros(Nz, Nx, Ny);
results_3D.T1MW = zeros(Nz, Nx, Ny);
results_3D.T1IEW = zeros(Nz, Nx, Ny);
results_3D.r2MW = zeros(Nz, Nx, Ny);
results_3D.r2IEW = zeros(Nz, Nx, Ny);
tic
for slc = 50:230
    y = img_cor(slc,:,:,:);
    mask_slc = mask_brain(slc,:,:);

    mask = mask_slc;

    % itinital values
    totalField_slc = totalField(slc,:,:);
    pini_slc = pini(slc,:,:);
    b1_slc = b1(slc,:,:);
    scaleFactor = max(PD_map(:));
    % PD_slc = PD_map(slc,:,:)/scaleFactor;
    PD_slc = PD_map(slc,:,:);
    PV1_slc = PV1_map(slc,:,:);
    PV2_slc = PV2_map(slc,:,:);

    pars0 = [];
    pars0.fieldmap   = totalField_slc;   % Hz    % it is not neccessary to provide the total field for initial starting point, yet it can provide the fitting results
    pars0.pini       = pini_slc;% rad
    pars0.b1map     = b1_slc;
    % [pars0.PDmap, ~] =  xy2voxels(PD_slc,mask_slc)/scaleFactor;
    pars0.fmw_init = PV1_slc;
    pars0.fiew_init = PV2_slc;
    

    if sum(mask(:))~=0


     objGPU                     = gpuMIMOSAMWI_A05(TE.',fixed_params);

    fitting = [];
    fitting.param = param;
    fitting.iteration           = 1000;
    fitting.initialLearnRate    = 0.01;% 0.01 ; 0.003
    fitting.decayRate           = 0;
    fitting.convergenceValue    = 1e-6;% slope
    fitting.tol                 = 1e-6;% noise level
    fitting.isDisplay           = 0;
    fitting.weightPower         = 0.5;
    % fitting.patience          = 50;     % allow longer plateau before stopping; sometimes fitting of all voxels are not finished yet
    % fitting.convergenceWindow = 40;     % larger window when smaller convengen => less sensitive to tiny oscillation
    % fitting.decayRate         = 1e-2;% to reduce the lr
    fitting.lambda = [3e-4];
    fitting.regmap = {'MWF'};
    fitting.TVmode = '3D';
    
    % yc, add
    fitting.isWeighted = false;
    fitting.debug         = 0;
    fitting.enableComplex = false;
    fitting.start = 'prior_mimosa';
    fitting.isComplex = 1;
    
    fitting.DIMWI.isFitIWF      = 1;
    fitting.DIMWI.isFitFreqMW   = 1;
    fitting.DIMWI.isFitFreqIEW   = 1;
    fitting.DIMWI.isFitR2sIEW    = 1;
    
    fitting.lossFunction = 'l1';   % deault L1 
    
    extraData = [];
    extraData.freqBKG   = pars0.fieldmap / (gpuMIMOSAMWI.gyro*fixed_params.B0); % in ppm
    extraData.pini      =  pars0.pini;
    extraData.IWF       = pars0.fiew_init;
    
    scale = prctile(PD_slc(:),98);

    extraData.PD       = PD_slc/scale;
    extraData.MWF       = PV1_slc;
    extraData.b1       = b1_slc;
    fitting.usingANN = true;   % deault L1 
    fitting.model = 'invivo';
    
    [out_ann]    = objGPU.estimate(y/scale, mask, extraData, fitting);


    results_3D.S0(slc, :, :) = sq(out_ann.final.S0);
    results_3D.MWF(slc, :, :) = sq(out_ann.final.MWF);
    results_3D.IEW(slc, :, :) = sq(out_ann.final.IEW);
    results_3D.FW(slc, :, :) = 1 - sq(out_ann.final.MWF + out_ann.final.IEW);
    
    results_3D.r2sMW(slc, :, :) = sq(out_ann.final.r2sMW);
    results_3D.r2sIEW(slc, :, :) = sq(out_ann.final.r2sIEW);
    
    results_3D.freqMW(slc, :, :) = sq(out_ann.final.freqMW * gpuMIMOSAMWI.gyro * fixed_params.B0);
    results_3D.freqIEW(slc, :, :) = sq(out_ann.final.freqIEW * gpuMIMOSAMWI.gyro * fixed_params.B0);
    results_3D.dfreqBKG(slc, :, :) = sq(out_ann.final.dfreqBKG * gpuMIMOSAMWI.gyro * fixed_params.B0);
    results_3D.dpini(slc, :, :) = sq(out_ann.final.dpini);
    
    results_3D.fieldmap_total(slc, :, :) = sq(out_ann.final.dfreqBKG * gpuMIMOSAMWI.gyro * fixed_params.B0 + pars0.fieldmap);
    results_3D.pini_total(slc, :, :) = sq(out_ann.final.dpini + pars0.pini);
    
    results_3D.T1MW(slc, :, :) = sq(out_ann.final.T1MW * 1000);
    results_3D.T1IEW(slc, :, :) = sq(out_ann.final.T1IEW * 1000);
    results_3D.r2MW(slc, :, :) = sq(out_ann.final.r2MW);
    results_3D.r2IEW(slc, :, :) = sq(out_ann.final.r2IEW);
    end
    
end
%%
save('gacelle_sub5_A05.mat','results_3D','-v7.3')
%%
MWF = single(sq( results_3D.MWF));
nii = make_nii(MWF);   
nii.hdr.dime.datatype = 16;          
nii.hdr.dime.bitpix = 32;
save_nii(nii,'mimosa_sub5_R8_MWF.nii')
%%
IEW = single(sq(results_3D.IEW));
nii = make_nii(IEW);   
nii.hdr.dime.datatype = 16;          
nii.hdr.dime.bitpix = 32;
save_nii(nii,'mimosa_sub5_R8_IEW.nii')
%%
FW = single(results_3D.FW);
nii = make_nii(FW);   
nii.hdr.dime.datatype = 16;         
nii.hdr.dime.bitpix = 32;
save_nii(nii,'mimosa_sub5_R8_FW.nii')
%%
%% ----------------------------chi sep-----------------------
%% χ-separation Tool

% This tool is MATLAB-based software forseparating para- and dia-magnetic susceptibility sources (χ-separation). 
% Separating paramagnetic (e.g., iron) and diamagnetic (e.g., myelin) susceptibility sources 
% co-existing in a voxel provides the distributions of two sources that QSM does not provides. 

% χ-separation tool v1.0

% Contact E-mail: snu.list.software@gmail.com 

% Reference
% H.-G. Shin, J. Lee, Y. H. Yun, S. H. Yoo, J. Jang, S.-H. Oh, Y. Nam, S. Jung, S. Kim, F. Masaki, W. 
% Kim, H. J. Choi, J. Lee. χ-separation: Magnetic susceptibility source separation toward iron and 
% myelin mapping in the brain. Neuroimage, 2021 Oct; 240:118371.

% χ-separation tool is powered by MEDI toolbox (for BET), STI Suite (for V-SHARP), SEGUE toolbox (for SEGUE), and mritools (for ROMEO).

%% Necessary preparation

% Set x-separation tool directory path
home_directory ='chi-separation_v2/chi-separation-main/Chisep_Toolbox_v1.1.2';
addpath(genpath(home_directory))

% Set MATLAB tool directory path 
% Set MATLAB tool directory path 
% xiangruili/dicm2nii (https://kr.mathworks.com/matlabcentral/fileexchange/42997-xiangruili-dicm2nii)

% Tools for NIfTI and ANALYZE image (https://kr.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)

% Download onnxconverter Add-on, and then install it.
% Deep Learning Toolbox Converter for ONNX Model Format 
% (https://kr.mathworks.com/matlabcentral/fileexchange/67296-deep-learning-toolbox-converter-for-onnx-model-format)

% Set QSM tool directory path 
% STI Suite (Version 3.0) (https://people.eecs.berkeley.edu/~chunlei.liu/software.html)

% MEDI toolbox (http://pre.weill.cornell.edu/mri/pages/qsm.html)

% SEGUE toolbox (https://xip.uclb.com/product/SEGUE)

% mritools toolbox (https://github.com/korbinian90/CompileMRI.jl/releases)

%% 
%% Run options - User define
RunOptions = struct();
% 'dicom': input DICOM | 'nifti': input NIfTI | Else: custom input (.mat)
RunOptions.InputType = 'Else';

% 'multi': multiple subjects | 'single: single-subject
RunOptions.multi = 'single';

% true: input brain mask | false: calculate brain mask
RunOptions.Mask = false;

% 'MEDI': MEDI brain extraction | 'custom': customize using FSL BET
RunOptions.Mask_method = 'MEDI';

% 'ARLO' | 'NNLS fitting' | 'Use preprocessed R2* or R2'' map'
RunOptions.R2sfit = 'ARLO'; 

% 'ROMEO + weighted echo averaging' | 'nonlinear complex fitting + SEGUE' | 'Laplacian'
RunOptions.Unwrap = 'ROMEO + weighted echo averaging'; 

% 'V-SHARP'
RunOptions.BFR = 'V-SHARP';

% 'Chi-sepnet' | 'Chi-separation (MEDI)' | 'Chi-separation (iLSQR)' 
RunOptions.Chisep = 'Chi-sepnet'; 

% 'Deep-learning' | 'Region-growing' | 'No'
RunOptions.VesselSeg = 'Deep-learning';

% GRE smoothing: 0 ~ 0.4(Default)
RunOptions.Tukey = double(0.4);

% 0: No inverse(Default) | 1: Inverse
RunOptions.PhaseInverse = 0;

% 1: have R2' | 0: don't have R2'
RunOptions.HaveR2Prime = 1;
% r2prime - R2' map in Hz unit (x, y, z). If you don't have R2' map, use chi-sepnet-R2* which doesn't require R2' map.

% 0: generate R2' from R2* using R2pnet | 1: generate R2' from R2* using scaling
RunOptions.is_scaling = 0;
RunOptions.scaling_factor = 0.19;

% false: No denoising for R2s | true: denosing for R2s
RunOptions.denoising = false;

% true: use resolution generalization | false: don't use
RunOptions.resgen = true; 
% Determine whether to use resolution generalization pipeline or to interpolate to 1 mm isotropic resolution
% 7T processing is available only with resolution generalization

RunOptions.OutputPath = '';
% Output path must not contatin ' '(spaces)

% Interpolation options (for B0 direction, Resampling)
% 'sinc' | 'spline'
RunOptions.interp_method = 'sinc';
RunOptions.sinc_window_size = 15;
% 'hann' | 'hamming' | 'blackman'
RunOptions.sinc_window_type = 'hann';

% Last stage Tukey
RunOptions.tukey_strength = 0.5;
RunOptions.tukey_pad = 0.1; %Recommend not to fix this

Data = struct();
Data.RunOptions = RunOptions;

meas = img_zsssl(:,:,:,end-9:end); clear img_zsssl;
meas = flip(permute(meas, [2 3 1 4]),1);

Data.CF = double(123138407);
Data.TE = 1.98:2.58:26.5
Data.B0dir = [0, 0, 1];
Data.VoxelSize = [1.0 1.0 1.0]
Data.Necho = size(meas,4);
Data.MatrixSize = size(meas);
Data.B0_strength = 3;
Data.MGRE_Mag = double(abs(meas));
Data.MGRE_Phs = double(angle(meas));
clearvars -except Params Data type_dir subj subj_dir path type type_path RunOptions home_directory

%% Fill in necessary parameters if empty
% Data.TE = [];                     % [ms]  [row vector]
% Data.B0dir = [];                  % []    [row vector]
% Data.CF = [];                     % [Hz]
% Data.B0_strength = [];            % [T]   B0_strength = CF / 42.58e6;


%% Params_check

% Force even dimension
input_field = {'MGRE_Mag','MGRE_Phs'};
for i = 1:length(input_field)
    [Data.(cell2mat(input_field(i))),x_odd,y_odd,z_odd] = even_pad(Data.(cell2mat(input_field(i))));
end
RunOptions.EvenSizePadding = [x_odd,y_odd,z_odd];
Data.MatrixSize = size(Data.MGRE_Mag);

% TE shape correction
if size(Data.TE,2) > 1
    Data.TE = Data.TE';
end

% Vendor options
if isfield(Data, 'Vendor') && strcmp(Data.Vendor,'P')
    RunOptions.Tukey = double(0);
    Data.RunOptions = RunOptions;
end
if isfield(Data, 'Vendor') && (strcmp(Data.Vendor,'G') || strcmp(Data.Vendor,'S'))
    Data.RunOptions = RunOptions;
end
    RunOptions.PhaseInverse = 0;



%% Tukey windowing
imgc = Data.MGRE_Mag .* exp(1i*Data.MGRE_Phs * (-1)^(RunOptions.PhaseInverse));
imgc = tukey_windowing(imgc,RunOptions.Tukey);
Data.MGRE_Mag_Tukey = abs(imgc);
Data.MGRE_Phs_Tukey = angle(imgc);

clearvars imgc


%% Brain mask (Range [0,1])
disp("=================< Brain masking >=================")
if RunOptions.Mask
    Data.Mask = load('mask.mat');
else
    if strcmp(RunOptions.Mask_method,'MEDI')                                % Use MEDI BET
        Data.Mask = BET(Data.MGRE_Mag_Tukey(:,:,:,1), Data.MatrixSize(1:3), Data.VoxelSize, 0.4);
%         Data.Mask = double(imerode(Data.Mask, strel('sphere',2)));
        Data.Mask = double(Data.Mask);
    else                                                                    % customizing using FSL
        mat2nii_ungz(Data.MGRE_Mag_Tukey,[Data.output_root,'\mag_tmp'])
        cmd = ['/home/user/fsl/bin/bet ',[Data.output_root,'\mag_tmp '],[Data.output_root,'\BET'], ' -m -R -f 0.55 -g 0.15 -S'];%-f 0.7 -g -0.08
        [status, result] = system(fsl_PathCorr(cmd));
        mask_brain = fliplr(rot90(niftiread([Data.output_root,'\BET_mask.nii.gz'])));
        Data.Mask = imerode(imdilate(mask_brain,strel('sphere',2)),strel('sphere',4));
    end
end

clearvars mask_brain
load('mapping_mwf_mimosa_R8_bp_zsssl.mat', 'T2_map', 'T2s_map')

R2s = 1./T2s_map*1000;
R2s(T2s_map==0)=0;
R2 = 1./T2_map*1000;
R2(T2_map==0)=0;
R2p = R2s - R2;
R2p(R2p<0)=0;

Data.R2p = flip(permute(R2p, [2 3 1]),1);
Data.R2p = Data.R2p(:,:,49:258);


if RunOptions.HaveR2Prime                                                   % Use Chi-sepnet-R2'
    Data.map = Data.R2p;
else                                                                        % Use Chi-sepnet-R2*
    Data.map = Data.R2s;
end
Data.map(Data.map < 0) = 0;


%% Calculate & correct bias field for Philips data
if isfield(Data,'Vendor')
    if(Data.Vendor == 'P')
        disp('Detecting bias field for Philips data')
        [biasField, detected] =  CustomBiasCorrection_step1(Data.MGRE_Phs_Tukey,logical(Data.Mask),Data.MGRE_Mag_Tukey);
        if detected
            Data.MGRE_Phs_BiasCor = CustomBiasCorrection_step2(Data.MGRE_Phs_Tukey,biasField);
            Data.RunOptions.PhilipsBiasCor = true;
        end
    end
end
if (isfield(Data,'MGRE_Phs_BiasCor'))
    phase = Data.MGRE_Phs_BiasCor;
elseif(isfield(Data,'MGRE_Phs_Tukey'))
    phase = Data.MGRE_Phs_Tukey;
else
    phase = Data.MGRE_Phs;
end


%% Phase Unwrapping (Range [-10,10] [rad])
% [unwrapped_phase[w*TE, angle]-> Echo combine -> UnwrappedPhase[w*dTE, angle]]
disp("================< Phase unwrapping >===============")
if(strcmp(RunOptions.Unwrap,'ROMEO + weighted echo averaging'))
    parameters.TE = Data.TE;
    parameters.mag = Data.MGRE_Mag_Tukey;
    parameters.mask = double(Data.Mask);
    parameters.calculate_B0 = false;
    parameters.phase_offset_correction = 'off';
    parameters.voxel_size = Data.VoxelSize;
    parameters.additional_flags = '-q -i';
    parameters.output_dir = ['romeo_tmp'];
    mkdir(parameters.output_dir);

    [unwrapped_phase, B0] = ROMEO(double(phase), parameters);
    unwrapped_phase(isnan(unwrapped_phase))= 0;

    % % Weighted echo averaging
    % TE_s = Data.TE/1000;
    % t2s_roi = 0.04;
    % W = (TE_s).*exp(-(TE_s)/t2s_roi);
    % weightedSum = 0;
    % TE_eff = 0;
    % for echo = 1:size(unwrapped_phase,4)
    %     weightedSum = weightedSum + W(echo)*unwrapped_phase(:,:,:,echo)./sum(W);
    %     TE_eff = TE_eff + W(echo)*TE_s(echo)./sum(W);
    % end
    % 
    % Data.UnwrappedPhase = weightedSum / TE_eff * (TE_s(2)-TE_s(1)) .* Data.Mask;

elseif(strcmp(RunOptions.Unwrap,'nonlinear complex fitting + SEGUE'))
    % Complex fitting from MEDI
    [field, error, residual_, phase0]=Fit_ppm_complex_TE(Data.MGRE_Mag_Tukey.*exp(-1i*phase), Data.TE);
   
    Inputs.Mask = double(Data.Mask); % 3D binary tissue mask, same size as one phase image
    Inputs.Phase = double(field); % For opposite phase
    Data.UnwrappedPhase = SEGUE(Inputs) .* Data.Mask; % Tissue phase in rad

elseif(strcmp(RunOptions.Unwrap,'Laplacian'))
    % Weighted echo combine + Laplacian
    [phase, N_std] = Preprocessing4Phase(Data.MGRE_Mag_Tukey,Data.MGRE_Phs_Tukey);
    pad_size=[12 12 12];
    [Data.UnwrappedPhase_, ~] = MRPhaseUnwrap(phase,'voxelsize',Data.VoxelSize,'padsize',pad_size);
    Data.UnwrappedPhase = Data.UnwrappedPhase / Data.dTE;

    % Laplacian + Echo sum
    pad_size=[12 12 12];
    [Data.UnwrappedPhase_, ~] = MRPhaseUnwrap(Data.MGRE_Phs_Tukey,'voxelsize',Data.VoxelSize,'padsize',pad_size);
    Data.UnwrappedPhase = sum(Data.UnwrappedPhase_,4) / sum(Data.TE);
   
    clearvars field_map pad_size
end
Data.local_field = 0;
iter_OE = 2;
for iter_ = 1:iter_OE
    te_gre_t = Data.TE(iter_:iter_OE:end)/1000;
    sc_f = Data.delta_TE/(te_gre_t(2) - te_gre_t(1));
    iFreq_sum = sum(unwrapped_phase(:, :, :, iter_:iter_OE:end), 4) / sum(te_gre_t) * (te_gre_t(2) - te_gre_t(1));
    [loc_f_sharp, Data.mask_brain_new] = V_SHARP(iFreq_sum, Data.Mask, 'voxelsize', Data.VoxelSize, 'smvsize', 20);
    Data.local_field = Data.local_field + loc_f_sharp/iter_OE*sc_f;
end
Data.delta_TE = (Data.TE(2)-Data.TE(1))/1000;
Data.local_field_hz = double(Data.local_field) / (2*pi*Data.delta_TE); % rad to hz
% %% Background field removal (Range [-5,5])
% % [local_field_hz [hz]]
% disp("============< Background field removal >============")
% if(strcmp(RunOptions.BFR,'V-SHARP'))
%     [Data.local_field, Data.mask_brain_new]=V_SHARP(Data.UnwrappedPhase, Data.Mask,'voxelsize', Data.VoxelSize,'smvsize', 20);
%     Data.delta_TE = (Data.TE(2)-Data.TE(1))/1000;
%     Data.local_field_hz = double(Data.local_field) / (2*pi*Data.delta_TE); % rad to hz
% end

%% QSM
% % 1. iLSQR from STI Suite
pad_size = [12, 12, 12];
Data.QSM = QSM_iLSQR(Data.local_field, Data.mask_brain_new,'TE',Data.delta_TE*1e3,'B0',Data.B0_strength,'H',Data.B0dir','padsize',pad_size,'voxelsize',Data.VoxelSize);

%% Chi separation
disp("============< χ-separation processing >============")
switch RunOptions.Chisep
    case 'Chi-sepnet'
        Dr = 114; % This parameter is different from the original paper (Dr = 137) because the network is trained on COSMOS-reconstructed maps
        if RunOptions.resgen
            % Use the resolution generalization pipeline. Resolution of input data is retained in the resulting chi-separation maps
            [Data.x_para, Data.x_dia, Data.x_tot, Data.qsm_map, Data.r2p_map] = chi_sepnet_general_new_wResolGen(home_directory, Data.local_field_hz, Data.map, Data.mask_brain_new, Dr, ...
                Data.B0dir, Data.CF, Data.VoxelSize, RunOptions.HaveR2Prime, Data.B0_strength, RunOptions.is_scaling, RunOptions.scaling_factor, RunOptions.interp_method, RunOptions.sinc_window_size, RunOptions.sinc_window_type);
        else
            % Interpolate the input maps to 1 mm isotropic resolution. 
            [Data.x_para, Data.x_dia, Data.x_tot, Data.qsm_map, Data.r2p_map] = chi_sepnet_general_sinc(home_directory, Data.local_field_hz, Data.map, Data.mask_brain_new, Dr, ...
                Data.B0dir, Data.CF, Data.VoxelSize, RunOptions.HaveR2Prime, Data.B0_strength, RunOptions.is_scaling, RunOptions.scaling_factor, RunOptions.interp_method, RunOptions.sinc_window_size, RunOptions.sinc_window_type);
        end
    
    case 'Chi-separation (MEDI)'
        Data.mag = sqrt(sum(Data.MGRE_Mag_Tukey.^2,4)) .* Data.mask_brain_new;
        Data.local_field_hz = Data.local_field_hz .* Data.mask_brain_new;
        Data.r2prime = Data.map .* Data.mask_brain_new;
        [~, N_std] = Preprocessing4Phase(Data.MGRE_Mag_Tukey, Data.MGRE_Phs_Tukey);
        params.b0_dir = Data.B0dir;
        params.CF = Data.CF;
        params.voxel_size = Data.VoxelSize;
        params.TE = Data.TE;
        params.lambda = 1;
        params.lambda_CSF = 1;
        params.Dr = 137;
        option_data.qsm = Data.QSM;
        option_data.mask_CSF = Data.mask_CSF;
        option_data.N_std = N_std;
        option_data.wG = [];
        option_data.wG_r2p = [];
        option_data.mask_FastRelax = zeros(size(Data.r2prime));
        option_data.mask_SlowRelax = zeros(size(Data.r2prime));
        [Data.x_para, Data.x_dia, Data.x_tot] = chi_sep_MEDI(Data.mag, Data.local_field_hz, Data.r2prime, N_std, Data.mask_brain_new, params, option_data);

    case 'Chi-separation (iLSQR)'
        Data.mag = sqrt(sum(Data.MGRE_Mag_Tukey.^2,4)) .* Data.mask_brain_new;
        Data.local_field_hz = Data.local_field_hz .* Data.mask_brain_new;
        Data.r2prime = Data.map .* Data.mask_brain_new;
        [~, N_std] = Preprocessing4Phase(Data.MGRE_Mag_Tukey, Data.MGRE_Phs_Tukey);
        params.b0_dir = Data.B0dir;
        params.CF = Data.CF;
        params.voxel_size = Data.VoxelSize;
        params.Dr = 137;
        option_data.qsm = Data.QSM;
        option_data.N_std = N_std;
        [Data.x_para, Data.x_dia, Data.x_tot] = chi_sep_iLSQR(Data.mag, Data.local_field_hz, Data.r2prime, Data.mask_brain_new, params, option_data);

end


if strcmp(RunOptions.interp_method, 'sinc')
    tukey_strength = RunOptions.tukey_strength;
    tukey_pad = RunOptions.tukey_pad;
    Data.x_para = real(tukey_windowing(Data.x_para,tukey_strength,round(size(Data.x_para).*tukey_pad))) .* Data.mask_brain_new;
    Data.x_dia = real(tukey_windowing(Data.x_dia,tukey_strength,round(size(Data.x_dia).*tukey_pad))) .* Data.mask_brain_new;
    Data.x_tot = real(tukey_windowing(Data.x_tot,tukey_strength,round(size(Data.x_tot).*tukey_pad))) .* Data.mask_brain_new;
    Data.qsm_map = real(tukey_windowing(Data.qsm_map,tukey_strength,round(size(Data.qsm_map).*tukey_pad))) .* Data.mask_brain_new;
    Data.r2p_map = real(tukey_windowing(Data.r2p_map,tukey_strength,round(size(Data.r2p_map).*tukey_pad))) .* Data.mask_brain_new;

    Data.x_para(Data.x_para < 0) = 0;
    Data.x_dia(Data.x_dia < 0) = 0;
    Data.r2p_map(Data.r2p_map < 0) = 0;
end

x_para = Data.x_para;
x_dia = Data.x_dia;
qsm = Data.qsm_map;
x_tot = Data.x_tot;

%%
save(['QSM_mwf_mimosa_R8_zsssl_sepnet.mat'],'x_para','x_dia','x_tot')
%%
if is_save_nii
nii_T1 = make_nii(sq((x_para)));
save_nii(nii_T1,'mimosa_R8_x_para.nii');
nii_T2 = make_nii(sq((x_dia)));
save_nii(nii_T2,'mimosa__R8_x_dia.nii');
nii_pd = make_nii(sq((x_tot)));
save_nii(nii_pd,'mimosa_R8_x_tot.nii');
end
