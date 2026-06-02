%% run_mimosa_feature_ablation.m
% English comments only

clear; close all; clc;

addpath(genpath('/autofs/cluster/berkin/yuting/MATLAB/demo/gacelle-main'));
addpath(genpath('/autofs/cluster/berkin/yuting/data/20251124_mwf_mimosa/MCR_mwf'));

%% ===================== Output settings =====================
root_prefix = 'MIMOSA_featAbl_loss1_tau0p01_noTemporal';
root_dir = fullfile(pwd, root_prefix);
if ~exist(root_dir, 'dir'); mkdir(root_dir); end

saveModels = true;   % set false to only save metrics/config
doPlotLoss = false;  % set true if you want per-run loss curve

%% ===================== Sequence parameters =====================
seq_param.esp             = 5.8e-3;
seq_param.turbo_factor    = 127;
seq_param.TR              = 4500e-3 - 5.8e-3*127 + 27.7e-3*127 - 162.4e-3*2;
seq_param.alpha_deg       = 4;
seq_param.num_reps        = 5;
seq_param.echo2use        = 1;
seq_param.gap_between_readouts  = 900e-3;
seq_param.time2relax_at_the_end = 0;
seq_param.ncontrast       = 14;

seq_param.TEs     = (1.98:2.58:26.5)*1e-3;
seq_param.TR_mte  = 27.7e-3;
seq_param.esp_mte = 2.58e-3;
seq_param.nechoes = numel(seq_param.TEs);

seq_param.TE_flash = 2.29e-3;
Npt = seq_param.ncontrast;

seq_param = build_MIMOSA14_acq_times(seq_param, seq_param.echo2use);

%% ===================== Fixed tissue constants =====================
T1_FW  = 4500e-3;
T2_FW  = 500e-3;
T2s_FW = 500e-3;
freq_FW_Hz = 0; 

%% ===================== Unknown ranges =====================
T1_MW_range_s   = [100e-3, 400e-3];
T1_IEW_range_s  = [700e-3, 1500e-3];

R2_MW_range_Hz  = [25, 100];
R2_IEW_range_Hz = [10, 17];

IEcoef_MW  = 0.5022;
IEcoef_IEW = 0.8567;
IEcoef_FW  = 0.9038;
% IEcoef_MW  = 1;
% IEcoef_IEW = 1;
% IEcoef_FW  = 1;

MWF_range       = [1e-8, 0.3];
IEW_range       = [1e-8, 1];

R2s_MW_cap_range_Hz  = [50, 200];
R2s_IEW_cap_range_Hz = [2,  50];

R2p_MW_range_Hz  = [ max(0, R2s_MW_cap_range_Hz(1)  - max(R2_MW_range_Hz)), ...
                     max(0, R2s_MW_cap_range_Hz(2)  - min(R2_MW_range_Hz)) ];

R2p_IEW_range_Hz = [ max(0, R2s_IEW_cap_range_Hz(1) - max(R2_IEW_range_Hz)), ...
                     max(0, R2s_IEW_cap_range_Hz(2) - min(R2_IEW_range_Hz)) ];

freq_MW_range_ppm   = [-0.05,  0.25];
freq_IEW_range_ppm  = [-0.10,  0.05];
freq_BKG_range_ppm  = [-0.40,  0.40];
ppm2Hz = 123.0489;

b1_range = [0.65, 1.35];

%% ===================== Network definition =====================
numOutputs = 2;
numLayers  = [];
numNeurons = [160 240 320 360 480 520 600]/2;

%% ===================== Training options =====================
executionEnvironment = "gpu";
Nepochs              = 1000;
alpha_leaky          = 0.01;

initialLearnRate = 0.01;
decayRate        = 0.0001;

cacheNSteps   = 20;
itersPerEpoch = 200;

%% ===================== Validation options =====================
Nvalid = 20000;

%% ===================== Build all ablation configs =====================
cfgs = build_feature_ablation_configs();
nCfg = numel(cfgs);

results = repmat(struct(), nCfg, 1);

%% ===================== Run all ablations =====================

for ic = 5%1:nCfg
    cfg = cfgs(ic);

    % English comment: reset RNG for fair comparison across configs
    rng(0, 'twister');

    % feature_idx = cfg.feature_idx(:).';
    feature_idx = [1:11];% wo time
    nfeatures   = numel(feature_idx);

    run_name = sprintf('A%02d_%s_F%02d', ic, cfg.name, nfeatures);
    out_dir  = fullfile(root_dir, run_name);
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end

    fprintf('\n==================== %s ====================\n', run_name);
    fprintf('Features: %s\n', mat2str(feature_idx));

    % English comment: create model + optimizer states
    parameters   = create_mlp(numLayers, numNeurons, nfeatures, numOutputs);
    averageGrad  = [];
    averageSqGrad = [];
    iteration = 0;

    % English comment: optional loss plot
    if doPlotLoss
        figLoss = figure('Name', ['Loss_' run_name]); %#ok<NASGU>
        C = colororder;
        lineLoss = animatedline('Color', C(2,:));
        grid on; xlabel('Iteration'); ylabel('Loss'); ylim([0 inf]); drawnow;
    end

    epoch_loss_mean = zeros(Nepochs, 1, 'single');

    % ===================== Training =====================
    for epoch = 1:Nepochs

        if epoch <= 40
            miniBatchSize = 128;
        elseif epoch <= 120
            miniBatchSize = 256;
        elseif epoch <= 300
            miniBatchSize = 512;
        else
            miniBatchSize = 1024;
        end

        % English comment: fill cache
        cache = cell(cacheNSteps,2);
        for c = 1:cacheNSteps
            [dlin3, dlout3] = make_minibatch_MIMOSA14_3pool( ...
                feature_idx, miniBatchSize, seq_param, ...
                T1_MW_range_s, R2_MW_range_Hz, ...
                T1_IEW_range_s, R2_IEW_range_Hz, ...
                T1_FW, T2_FW, T2s_FW, ...
                IEcoef_MW, IEcoef_IEW, IEcoef_FW, ...
                MWF_range, IEW_range, ...
                R2p_MW_range_Hz, R2p_IEW_range_Hz, ...
                R2s_MW_cap_range_Hz, R2s_IEW_cap_range_Hz, ...
                freq_MW_range_ppm, freq_IEW_range_ppm, freq_BKG_range_ppm, ppm2Hz, ...
                b1_range);

            if executionEnvironment == "gpu"
                dlin3  = gpuArray(dlin3);
                dlout3 = gpuArray(dlout3);
            end
            cache{c,1} = dlin3;
            cache{c,2} = dlout3;
        end
        cacheIdx = 1;

        epoch_loss_total = 0;

        for it = 1:itersPerEpoch
            iteration = iteration + 1;

            dlin3  = cache{cacheIdx,1};
            dlout3 = cache{cacheIdx,2};
            cacheIdx = cacheIdx + 1;
            if cacheIdx > cacheNSteps
                cacheIdx = 1;
            end

            dlin  = reshape(dlin3,  [nfeatures,  miniBatchSize*Npt]);
            dlout = reshape(dlout3, [numOutputs, miniBatchSize*Npt]);
            dlin  = dlarray(dlin,'CB');
            dlout = dlarray(dlout,'CB');

            [gradients, loss, loss_norm] = dlfeval( ...
                @modelGradients_mlp_MIMOSA_normOnly, parameters, dlin, dlout, alpha_leaky); %#ok<NASGU>

            learningRate = initialLearnRate / (1 + decayRate*iteration);
            [parameters, averageGrad, averageSqGrad] = adamupdate( ...
                parameters, gradients, averageGrad, averageSqGrad, iteration, learningRate);

            lossVal = double(gather(extractdata(loss)));
            epoch_loss_total = epoch_loss_total + lossVal;

            if doPlotLoss
                addpoints(lineLoss, iteration, lossVal);
                drawnow limitrate;
            end
        end

        epoch_loss_mean(epoch) = single(epoch_loss_total / itersPerEpoch);

        if mod(epoch,10) == 0
            fprintf('Run %s | Epoch %d | loss=%.6g | batch=%d\n', ...
                run_name, epoch, epoch_loss_mean(epoch), miniBatchSize);
        end
    end

    % English comment: pack trained net
    dlnet = struct();
    dlnet.parameters = parameters;
    dlnet.epoch      = Nepochs;
    dlnet.alpha      = alpha_leaky;

    if saveModels
        save(fullfile(out_dir, [run_name '_model.mat']), 'dlnet', 'epoch_loss_mean', 'feature_idx');
    else
        save(fullfile(out_dir, [run_name '_trainStats.mat']), 'epoch_loss_mean', 'feature_idx');
    end

    % ===================== Validation =====================
    [valMetrics, valPerPoint] = validate_model( ...
        dlnet.parameters, dlnet.alpha, feature_idx, executionEnvironment, ...
        Nvalid, seq_param, ...
        T1_MW_range_s, R2_MW_range_Hz, ...
        T1_IEW_range_s, R2_IEW_range_Hz, ...
        T1_FW, T2_FW, T2s_FW, ...
        IEcoef_MW, IEcoef_IEW, IEcoef_FW, ...
        MWF_range, IEW_range, ...
        R2p_MW_range_Hz, R2p_IEW_range_Hz, ...
        R2s_MW_cap_range_Hz, R2s_IEW_cap_range_Hz, ...
        freq_MW_range_ppm, freq_IEW_range_ppm, freq_BKG_range_ppm, ppm2Hz, ...
        b1_range);

    save(fullfile(out_dir, [run_name '_valMetrics.mat']), 'valMetrics', 'valPerPoint', 'feature_idx', 'cfg');

    % English comment: record summary
    results(ic).run_name       = run_name;
    results(ic).cfg_name       = cfg.name;
    results(ic).feature_idx    = feature_idx;
    results(ic).nfeatures      = nfeatures;
    results(ic).useFW          = cfg.useFW;
    results(ic).freq_mode      = cfg.freq_mode;
    results(ic).time_mode      = cfg.time_mode;
    results(ic).final_train_loss = double(epoch_loss_mean(end));
    results(ic).val_complex_nrse = valMetrics.complex_nrse;
    results(ic).val_complex_mae  = valMetrics.complex_mae;
    results(ic).val_mean_abs_dphi = valMetrics.mean_abs_dphi;
    results(ic).val_mag_mae      = valMetrics.mag_mae;

    % English comment: print key metrics
    fprintf('Run %s | Val Complex-NRSE=%.6g | MAE=%.6g | mean|dphi|=%.6g | MagMAE=%.6g\n', ...
        run_name, valMetrics.complex_nrse, valMetrics.complex_mae, valMetrics.mean_abs_dphi, valMetrics.mag_mae);
end

%% ===================== Save master results =====================
save(fullfile(root_dir, 'ALL_results.mat'), 'results', 'cfgs');

% English comment: also write a compact CSV summary
T = struct2table(results);
writetable(T, fullfile(root_dir, 'ALL_results_summary.csv'));

disp('All ablation runs finished and results saved.');

%% ========================================================================
%% Local functions
%% ========================================================================

function cfgs = build_feature_ablation_configs()
% English comments only

% Full feature channels (fixed mapping in generator):
% 1:MWF 2:IEW 3:T2sMWdecay 4:T2sIEWdecay 5:T1MWn 6:T1IEWn 7:R2MWn 8:R2IEWn 9:B1
% 10:FreqMW 11:FreqIEW
% 12:sinMW 13:cosMW 14:sinIEW 15:cosIEW
% 16:PT 17-20:prepOH 21:t_acq_norm 22:FW

mandatory = 1:9;

freqA = [10 11];           % raw freq
freqB = 12:15;             % sin/cos
freqSets = {freqA, freqB, [freqA freqB]};
freqNames = {'freqRaw', 'phaseTrig', 'freqRaw+phaseTrig'};

timeGroups = { ...
    [16], ...              % PT
    17:20, ...             % prepOH
    [21], ...               % t_acq_norm
    [],...
};
timeNames = {'PT', 'prepOH', 't_acq','None'};

% All non-empty subsets of the 3 time groups
timeSets = {};
timeSetNames = {};
for m = 1:7
    idx = find(bitget(m,1:3));
    tmp = [];
    nm  = {};
    for k = 1:numel(idx)
        tmp = [tmp timeGroups{idx(k)}]; %#ok<AGROW>
        nm{end+1} = timeNames{idx(k)}; %#ok<AGROW>
    end
    timeSets{end+1} = tmp; %#ok<AGROW>
    timeSetNames{end+1} = strjoin(nm, '+'); %#ok<AGROW>
end

cfgs = struct('name', {}, 'feature_idx', {}, 'freq_mode', {}, 'time_mode', {}, 'useFW', {});

ic = 0;
for ifr = 1:numel(freqSets)
    for it = 1:numel(timeSets)
        for useFW = 0:1
            ic = ic + 1;

            feat = [mandatory freqSets{ifr} timeSets{it}];
            if useFW == 1
                feat = [feat 22];
            end
            feat = unique(feat, 'stable');

            cfgs(ic).freq_mode = freqNames{ifr}; %#ok<AGROW>
            cfgs(ic).time_mode = timeSetNames{it}; %#ok<AGROW>
            cfgs(ic).useFW     = logical(useFW); %#ok<AGROW>
            cfgs(ic).feature_idx = feat; %#ok<AGROW>

            cfgs(ic).name = sprintf('%s__%s__FW%d', cfgs(ic).freq_mode, cfgs(ic).time_mode, useFW); %#ok<AGROW>
        end
    end
end
end

function [gradients, loss, loss_norm] = modelGradients_mlp_MIMOSA_normOnly(parameters, dlXf, dlRf, alpha)
% English comments only

epsA = dlarray(single(1e-4),'CB');
tau  = dlarray(single(1e-2),'CB');

if nargin < 4 || isempty(alpha)
    alpha = 0.01;
end

U = mlp_model_leakyRelu(parameters, dlXf, alpha); % [2,N]

Ur = dlRf(1,:);  Ui = dlRf(2,:);
Pr = U(1,:);     Pi = U(2,:);

mag_gt = sqrt(Ur.^2 + Ui.^2 + epsA);

wA = (mag_gt.^2) ./ (mag_gt.^2 + tau.^2);

dPr = Pr - Ur;
dPi = Pi - Ui;

err_c = sqrt(dPr.^2 + dPi.^2 + epsA);

loss_norm = mean(wA .* (err_c ./ (mag_gt + epsA)), "all");
loss = loss_norm;

gradients = dlgradient(loss, parameters);
end

function [valMetrics, perPoint] = validate_model( ...
    parameters, alpha_leaky, feature_idx, executionEnvironment, ...
    Nvalid, seq_param, ...
    T1_MW_range_s, R2_MW_range_Hz, ...
    T1_IEW_range_s, R2_IEW_range_Hz, ...
    T1_FW, T2_FW, T2s_FW, ...
    IEcoef_MW, IEcoef_IEW, IEcoef_FW, ...
    MWF_range, IEW_range, ...
    R2p_MW_range_Hz, R2p_IEW_range_Hz, ...
    R2s_MW_cap_range_Hz, R2s_IEW_cap_range_Hz, ...
    freq_MW_range_ppm, freq_IEW_range_ppm, freq_BKG_range_ppm, ppm2Hz, ...
    b1_range)

% English comments only

Npt = seq_param.ncontrast;
numOutputs = 2;

[dlin3, dlout3] = make_minibatch_MIMOSA14_3pool( ...
    feature_idx, Nvalid, seq_param, ...
    T1_MW_range_s, R2_MW_range_Hz, ...
    T1_IEW_range_s, R2_IEW_range_Hz, ...
    T1_FW, T2_FW, T2s_FW, ...
    IEcoef_MW, IEcoef_IEW, IEcoef_FW, ...
    MWF_range, IEW_range, ...
    R2p_MW_range_Hz, R2p_IEW_range_Hz, ...
    R2s_MW_cap_range_Hz, R2s_IEW_cap_range_Hz, ...
    freq_MW_range_ppm, freq_IEW_range_ppm, freq_BKG_range_ppm, ppm2Hz, ...
    b1_range);

nfeatures = numel(feature_idx);

dlin  = reshape(dlin3,  [nfeatures,  Nvalid*Npt]);
dlout = reshape(dlout3, [numOutputs, Nvalid*Npt]);

dlX = dlarray(dlin, 'CB');
dlR = dlarray(dlout,'CB');

if executionEnvironment == "gpu"
    dlX = gpuArray(dlX);
    dlR = gpuArray(dlR);
end

U = mlp_model_leakyRelu(parameters, dlX, alpha_leaky);
U = gather(extractdata(U));
R = gather(extractdata(dlR));

pred_reim = reshape(U, [2, Nvalid, Npt]);
gt_reim   = reshape(R, [2, Nvalid, Npt]);

pred_c = squeeze(pred_reim(1,:,:) + 1i*pred_reim(2,:,:));
gt_c   = squeeze(gt_reim(1,:,:)   + 1i*gt_reim(2,:,:));

err_abs = abs(pred_c - gt_c);
gt_abs  = abs(gt_c);

valMetrics = struct();

valMetrics.mean_abs_gt   = mean(gt_abs(:));
valMetrics.mean_abs_pred = mean(abs(pred_c(:)));
valMetrics.mean_abs_err  = mean(err_abs(:));

tmp = pred_c(:) .* conj(gt_c(:));
dphi = angle(tmp);

valMetrics.mean_dphi      = mean(dphi);
valMetrics.std_dphi       = std(dphi);
valMetrics.mean_abs_dphi  = mean(abs(dphi));

tmpErr = err_abs(:);
tmpGt  = gt_abs(:);
valMetrics.complex_nrse = sqrt(sum(tmpErr.^2)) / sqrt(max(sum(tmpGt.^2), eps('single')));
valMetrics.complex_mae  = mean(tmpErr);

re_err = real(pred_c(:)) - real(gt_c(:));
im_err = imag(pred_c(:)) - imag(gt_c(:));
re_gt  = real(gt_c(:));
im_gt  = imag(gt_c(:));

valMetrics.nrse_re = sqrt(sum(re_err.^2)) / sqrt(max(sum(re_gt.^2), eps('single')));
valMetrics.nrse_im = sqrt(sum(im_err.^2)) / sqrt(max(sum(im_gt.^2), eps('single')));

mag_err = abs(pred_c(:)) - abs(gt_c(:));
mag_gt  = abs(gt_c(:));
valMetrics.mag_nrse = sqrt(sum(mag_err.^2)) / sqrt(max(sum(mag_gt.^2), eps('single')));
valMetrics.mag_mae  = mean(abs(mag_err));

% English comment: per-point complex NRSE
nrse_pt = zeros(1, Npt, 'single');
for p = 1:Npt
    e = pred_c(:,p) - gt_c(:,p);
    e2 = sum(abs(e).^2);
    g2 = sum(abs(gt_c(:,p)).^2);
    nrse_pt(p) = sqrt(e2 / max(g2, eps('single')));
end

perPoint = struct();
perPoint.complex_nrse_pt = nrse_pt;
end

function [dlin3, dlout3] = make_minibatch_MIMOSA14_3pool(feature_idx, batchSize, seq, ...
    T1_MW_range_s, R2_MW_range_Hz, ...
    T1_IEW_range_s, R2_IEW_range_Hz, ...
    T1_FW, T2_FW, T2s_FW, ...
    IEcoef_MW, IEcoef_IEW, IEcoef_FW, ...
    MWF_range, IEW_range, ...
    R2p_MW_range_Hz, R2p_IEW_range_Hz, ...
    R2s_MW_cap_range_Hz, R2s_IEW_cap_range_Hz, ...
    freq_MW_range_ppm, freq_IEW_range_ppm, freq_BKG_range_ppm, ppm2Hz, ...
    b1_range)
% English comments only

Npt = seq.ncontrast;
numOutputs = 2;

% English comment: always build the full 22-channel feature tensor, then subset
nFull = 22;
featuresFull = zeros(batchSize, Npt, nFull, 'single');
label        = zeros(batchSize, Npt, numOutputs, 'single');

% ---- Sample fractions ----
MWF = single(rand(1,batchSize) * diff(MWF_range) + min(MWF_range));

IEW = zeros(1,batchSize,'single');
IEW_min = IEW_range(1);
for k = 1:batchSize
    IEW_max_k = min(IEW_range(2), 1 - MWF(k));
    if IEW_max_k < IEW_min
        IEW_max_k = IEW_min;
    end
    IEW(k) = single(rand() * (IEW_max_k - IEW_min) + IEW_min);
end

% ---- Sample T1/R2 ----
T1_MW_s  = single(rand(1,batchSize) * diff(T1_MW_range_s)  + min(T1_MW_range_s));
T1_IEW_s = single(rand(1,batchSize) * diff(T1_IEW_range_s) + min(T1_IEW_range_s));

R2_MW_Hz  = single(rand(1,batchSize) * diff(R2_MW_range_Hz)  + min(R2_MW_range_Hz));
R2_IEW_Hz = single(rand(1,batchSize) * diff(R2_IEW_range_Hz) + min(R2_IEW_range_Hz));

T2_MW_s  = 1 ./ (R2_MW_Hz  + eps('single'));
T2_IEW_s = 1 ./ (R2_IEW_Hz + eps('single'));

% ---- Sample R2p and derive R2s ----
R2pMW_min = max(0, R2s_MW_cap_range_Hz(1) - R2_MW_Hz);
R2pMW_max = max(0, R2s_MW_cap_range_Hz(2) - R2_MW_Hz);

R2pMW_min = max(R2pMW_min, R2p_MW_range_Hz(1));
R2pMW_max = min(R2pMW_max, R2p_MW_range_Hz(2));

badMW = (R2pMW_max < R2pMW_min);
R2pMW_max(badMW) = R2pMW_min(badMW);

uMW   = rand(1, batchSize, 'single');
R2pMW = single(R2pMW_min) + uMW .* single(R2pMW_max - R2pMW_min);

R2pIEW_min = max(0, R2s_IEW_cap_range_Hz(1) - R2_IEW_Hz);
R2pIEW_max = max(0, R2s_IEW_cap_range_Hz(2) - R2_IEW_Hz);

R2pIEW_min = max(R2pIEW_min, R2p_IEW_range_Hz(1));
R2pIEW_max = min(R2pIEW_max, R2p_IEW_range_Hz(2));

badIEW = (R2pIEW_max < R2pIEW_min);
R2pIEW_max(badIEW) = R2pIEW_min(badIEW);

uIEW   = rand(1, batchSize, 'single');
R2pIEW = single(R2pIEW_min) + uIEW .* single(R2pIEW_max - R2pIEW_min);

R2sMW  = R2_MW_Hz  + R2pMW;
R2sIEW = R2_IEW_Hz + R2pIEW;

% ---- Frequency and B1 ----
freqMW_ppm  = single(rand(1,batchSize) * diff(freq_MW_range_ppm)  + min(freq_MW_range_ppm));
freqIEW_ppm = single(rand(1,batchSize) * diff(freq_IEW_range_ppm) + min(freq_IEW_range_ppm));
freqBKG_ppm = single(rand(1,batchSize) * diff(freq_BKG_range_ppm) + min(freq_BKG_range_ppm));

b1 = single(rand(1,batchSize) * diff(b1_range) + min(b1_range));

% ---- Fractions ----
f_mw  = MWF;
f_iew = IEW;
f_fw  = 1 - MWF - IEW;

% ---- Local frequency (dffw0) ----
lf_mw_Hz  = freqMW_ppm  * ppm2Hz;
lf_iew_Hz = freqIEW_ppm * ppm2Hz;
lf_bkg_Hz = freqBKG_ppm * ppm2Hz; %#ok<NASGU>

local_mw_Hz  = lf_mw_Hz;
local_iew_Hz = lf_iew_Hz;
local_fw_Hz  = zeros(size(lf_iew_Hz));

% ---- TE list ----
TE_all = [seq.TE_flash seq.TE_flash seq.TE_flash seq.TE_flash seq.TEs(:)'];

% ---- Decay and phase ----
T2s_MW_decay  = exp(- (R2sMW(:)  * TE_all));
T2s_IEW_decay = exp(- (R2sIEW(:) * TE_all));

phi_mw_tot  = 2*pi * (local_mw_Hz(:)  * TE_all);
phi_iew_tot = 2*pi * (local_iew_Hz(:) * TE_all);

sinPhi_MW  = sin(phi_mw_tot);
cosPhi_MW  = cos(phi_mw_tot);
sinPhi_IEW = sin(phi_iew_tot);
cosPhi_IEW = cos(phi_iew_tot);

% ---- Simulate pools ----
t2s_mw  = 1 ./ (R2sMW(:)  + eps('single'));
t2s_iew = 1 ./ (R2sIEW(:) + eps('single'));
t2s_fw  = T2s_FW;

[~, Mxy_mtx_mw, ~] = sim_mwf_MIMOSA_bp_mxy0_cx_v2_dlarray( ...
    seq.TR, seq.alpha_deg, seq.esp, seq.turbo_factor, ...
    single(T1_MW_s(:)), single(T2_MW_s(:)), seq.num_reps, seq.echo2use, ...
    seq.TR_mte, seq.esp_mte, seq.TEs, ...
    single(t2s_mw(:)), seq.gap_between_readouts, seq.time2relax_at_the_end, ...
    b1(:), IEcoef_MW, local_mw_Hz(:), f_mw(:));

[~, Mxy_mtx_iew, ~] = sim_mwf_MIMOSA_bp_mxy0_cx_v2_dlarray( ...
    seq.TR, seq.alpha_deg, seq.esp, seq.turbo_factor, ...
    single(T1_IEW_s(:)), single(T2_IEW_s(:)), seq.num_reps, seq.echo2use, ...
    seq.TR_mte, seq.esp_mte, seq.TEs, ...
    single(t2s_iew(:)), seq.gap_between_readouts, seq.time2relax_at_the_end, ...
    b1(:), IEcoef_IEW, local_iew_Hz(:), f_iew(:));

[~, Mxy_mtx_fw, ~] = sim_mwf_MIMOSA_bp_mxy0_cx_v2_dlarray( ...
    seq.TR, seq.alpha_deg, seq.esp, seq.turbo_factor, ...
    single(T1_FW), single(T2_FW), seq.num_reps, seq.echo2use, ...
    seq.TR_mte, seq.esp_mte, seq.TEs, ...
    single(t2s_fw), seq.gap_between_readouts, seq.time2relax_at_the_end, ...
    b1(:), IEcoef_FW, local_fw_Hz(:), f_fw(:));

smw  = squeeze(Mxy_mtx_mw(:,:,end));
siew = squeeze(Mxy_mtx_iew(:,:,end));
sfw  = squeeze(Mxy_mtx_fw(:,:,end));

s_mix = smw + siew + sfw;

label(:,:,1) = single(real(s_mix)).';
label(:,:,2) = single(imag(s_mix)).';

% ---- Build feature tensors (full 22) ----
ptIndexNorm = single((0:Npt-1) / (Npt-1));

MWF_2d      = repmat(MWF(:), 1, Npt);
IEW_2d      = repmat(IEW(:), 1, Npt);
FW_2d       = repmat(f_fw(:), 1, Npt);

Freq_MW_2d  = repmat(freqMW_ppm(:), 1, Npt);
Freq_IEW_2d = repmat(freqIEW_ppm(:), 1, Npt);

B1_2d       = repmat(b1(:), 1, Npt);
PT_2d       = repmat(ptIndexNorm, batchSize, 1);

T1_MW_norm  = (T1_MW_s  - min(T1_MW_range_s))  ./ diff(T1_MW_range_s);
T1_IEW_norm = (T1_IEW_s - min(T1_IEW_range_s)) ./ diff(T1_IEW_range_s);
R2_MW_norm  = (R2_MW_Hz  - min(R2_MW_range_Hz))  ./ diff(R2_MW_range_Hz);
R2_IEW_norm = (R2_IEW_Hz - min(R2_IEW_range_Hz)) ./ diff(R2_IEW_range_Hz);

T1MWn_2d  = repmat(T1_MW_norm(:),  1, Npt);
T1IEWn_2d = repmat(T1_IEW_norm(:), 1, Npt);
R2MWn_2d  = repmat(R2_MW_norm(:),  1, Npt);
R2IEWn_2d = repmat(R2_IEW_norm(:), 1, Npt);

t_acq14_s = single(seq.t_acq14_s(:).');
t_acq_2d  = repmat(t_acq14_s, batchSize, 1);
tmax      = max(t_acq14_s);
t_acq_norm = t_acq_2d ./ (tmax + eps('single'));

% 1..9 mandatory
featuresFull(:,:,1) = MWF_2d;
featuresFull(:,:,2) = IEW_2d;
featuresFull(:,:,3) = single(T2s_MW_decay);
featuresFull(:,:,4) = single(T2s_IEW_decay);
featuresFull(:,:,5) = T1MWn_2d;
featuresFull(:,:,6) = T1IEWn_2d;
featuresFull(:,:,7) = R2MWn_2d;
featuresFull(:,:,8) = R2IEWn_2d;
featuresFull(:,:,9) = B1_2d;

% frequency group
featuresFull(:,:,10) = Freq_MW_2d;
featuresFull(:,:,11) = Freq_IEW_2d;

% phase trig group
featuresFull(:,:,12) = single(sinPhi_MW);
featuresFull(:,:,13) = single(cosPhi_MW);
featuresFull(:,:,14) = single(sinPhi_IEW);
featuresFull(:,:,15) = single(cosPhi_IEW);

% time group
featuresFull(:,:,16) = PT_2d;

prepOH = zeros(batchSize, Npt, 4, 'like', featuresFull);
prepOH(:,1,1)=1; prepOH(:,2,2)=1; prepOH(:,3,3)=1; prepOH(:,4,4)=1;
featuresFull(:,:,17:20) = prepOH;

featuresFull(:,:,21) = t_acq_norm;

% FW optional
featuresFull(:,:,22) = FW_2d;

% ---- Subset ----
features = featuresFull(:,:,feature_idx);

dlin3  = permute(features, [3 1 2]);
dlout3 = permute(label,    [3 1 2]);
end