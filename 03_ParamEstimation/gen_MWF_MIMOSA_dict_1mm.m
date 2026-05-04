%% generate dict, only once 
%% B1 
thre_high   = 1.35;
thre_low    = 0.65;
num_b1_bins = 50; % default: 50  --> decrease to 20
b1_val      = linspace( thre_low, thre_high, num_b1_bins );
%--------------------------------------------------------------------------
%% create look up table 
%--------------------------------------------------------------------------
tic
clc
fprintf('generating dictionaries...\n');

t1_entries  = [5:10:3000, 3100:50:5000];
t2_entries  = [1:2:350, 370:20:1000, 1100:100:3000];
t2s_entries  = [1:1:100, 105:5:200 210:50:500];


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

t1t2_lut_prune = zeros( [length(t1t2_lut) - idx, 2] );

%% load look up table
addpath(genpath('utils/'))
load('dict/ielookup_4mimosa.mat');
ie_lkp = reshape((ielookup.ies_mtx),[length(t1t2_lut),1]);
ie_lkp_prune = zeros( [length(t1t2_lut) - idx, 1] );
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
%--------------------------------------------------------------------------
%% set qalas 4acq mgre d1 param 
%--------------------------------------------------------------------------
param.esp             = 5.8 * 1e-3;
param.turbo_factor    = 127;
param.TR          = 4500e-3 - 5.8*127e-3 + 27.7*127e-3 - 162.4*2e-3;
param.alpha_deg   = 4;
param.num_reps    = 5;
param.echo2use    = 1;
param.gap_between_readouts    = 900e-3;
param.time2relax_at_the_end   = 0;

param.nconsrast = 14;

% ##### mgre
param.TE_mte = [1.98:2.58:26.5].*1e-3;
param.TR_mte = 27.7e-3;
param.esp_mte = 2.58e-3;
nechoes = length(param.TE_mte);


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

%% 
T2s_val = t2s_entries;


signal = zeros([length(t1t2_lut_prune),14, length(T2s_val),  length(b1_val)]);

length_t2s_val   = length(T2s_val);
length_b1_val   = length(b1_val);

cnt             = 0;
iniTime         = clock;
iniTime0        = iniTime;

% % parallel computing
delete(gcp('nocreate'))
c = parcluster('local');
total_cores = c.NumWorkers;
parpool(4)

parfor t2s = 1:length_t2s_val
    for b1 = 1:length_b1_val
            cnt = cnt + 1;
            
            [Mz, Mxy] = sim_mwf_MIMOSA_bp(TR, alpha_deg, esp, turbo_factor, t1t2_lut_prune(:,1)*1e-3, t1t2_lut_prune(:,2)*1e-3, num_reps, echo2use,TR_mte,esp_mte,TEs, T2s_val(t2s)*1e-3, gap_between_readouts, time2relax_at_the_end, b1_val(b1), ie_lkp_prune(:,1));
            
            temp = abs(Mxy(:,:,end).');
            
            for n = 1:size(temp,1)
                temp(n,:) = temp(n,:) / sum(abs(temp(n,:)).^2)^0.5;
            end
            
            signal(:,:,t2s,b1) = temp;
            
    end
end
% delete(gcp('nocreate'))
fprintf('total elapsed time: %.1f sec\n\n',etime(clock,iniTime0));


length_dict = length(t1t2_lut_prune);

dict = signal;
clear signal
% save('dict_mwf_mimosa_bp_50b1_126T2s.mat','dict','-v7.3')
%%
h5create('dict/dict_mwf_mimosa_bp_50b1_126T2s_double.h5', '/dict', size(dict), ...
         'Datatype', 'double', ...
         'Chunksize', [1024 14 10 1], ... 
         'Deflate', 6);
h5write('dict/dict_mwf_mimosa_bp_50b1_126T2s_double.h5', '/dict', dict);