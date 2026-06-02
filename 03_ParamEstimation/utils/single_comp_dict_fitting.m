T1_all = zeros(size(img_all,1),size(img_all,2),size(img_all,3));
T2_all = zeros(size(img_all,1),size(img_all,2),size(img_all,3));
PD_all = zeros(size(img_all,1),size(img_all,2),size(img_all,3));
T2s_all = zeros(size(img_all,1),size(img_all,2),size(img_all,3));
IE_all = zeros(size(img_all,1),size(img_all,2),size(img_all,3));

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

% figure;
% plot(values, counts);

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
