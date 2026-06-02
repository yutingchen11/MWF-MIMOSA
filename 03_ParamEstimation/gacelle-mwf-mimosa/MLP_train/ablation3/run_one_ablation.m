function run_one_ablation(cfg_name, seed)
% English comments only
%
% Single-ablation training + validation entry point (arch + feature ablation).
% - numNeurons: base / {1,2,4,8}
% - features: 1-9 mandatory + (freqRaw OR phaseTrig OR both) + fixed PT+prepOH
% - NO FW feature (22)
%
% Usage:
%   run_one_ablation('NDIV2__freqRaw__PT+prepOH')
%   run_one_ablation(1)
%   run_one_ablation('NDIV4__phaseTrig__PT+prepOH', 0)

%% ---------------- Input parsing ----------------
if nargin < 1 || isempty(cfg_name)
    error('cfg_name is required.');
end
if nargin < 2 || isempty(seed)
    seed = 0;
end

% English comment: deterministic RNG for reproducibility
rng(seed, 'twister');
try
    parallel.gpu.rng(seed, 'Philox4x32-10');
catch
end

cfgs = build_ablation_configs_arch_feat();
cfg_id = resolve_cfg_id(cfg_name, cfgs);
cfg = cfgs(cfg_id);

feature_idx = cfg.feature_idx(:).';
nfeatures   = numel(feature_idx);

%% ---------------- Paths ----------------
% English comment: addpath is not allowed in deployed (compiled) mode
if ~isdeployed
    addpath(genpath('/autofs/cluster/berkin/yuting/MATLAB/demo/gacelle-main'));
    addpath(genpath('/autofs/cluster/berkin/yuting/data/20251124_mwf_mimosa/MCR_mwf'));
end

root_prefix = 'MIMOSA_featAbl_ARCHxFEAT_loss1_tau0p01_v1';
root_dir = fullfile(pwd, root_prefix);
if ~exist(root_dir, 'dir')
    mkdir(root_dir);
end

run_name = sprintf('A%02d_%s_F%02d_seed%d', cfg_id, sanitize_name(cfg.name), nfeatures, seed);
out_dir  = fullfile(root_dir, run_name);
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

save(fullfile(out_dir, [run_name '_config.mat']), 'cfg', 'cfg_id', 'feature_idx', 'seed');

fprintf('\n============================================================\n');
fprintf('Running ablation: %s\n', cfg.name);
fprintf('Config ID       : %d\n', cfg_id);
fprintf('Seed            : %d\n', seed);
fprintf('Feature count   : %d\n', nfeatures);
fprintf('Feature idx     : %s\n', mat2str(feature_idx));
fprintf('Neuron div      : %d\n', cfg.neuron_div);
fprintf('numNeurons      : %s\n', mat2str(cfg.numNeurons));
fprintf('Output dir      : %s\n', out_dir);
fprintf('============================================================\n');

%% ---------------- Sequence parameters ----------------
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

%% ---------------- Fixed tissue constants ----------------
T1_FW  = 4500e-3;
T2_FW  = 500e-3;
T2s_FW = 500e-3;

%% ---------------- Unknown ranges ----------------
T1_MW_range_s   = [100e-3, 400e-3];
T1_IEW_range_s  = [700e-3, 1500e-3];

R2_MW_range_Hz  = [25, 100];
R2_IEW_range_Hz = [10, 17];

IEcoef_MW  = 0.5022;
IEcoef_IEW = 0.8567;
IEcoef_FW  = 0.9038;

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

%% ---------------- Network definition ----------------
numOutputs = 2;
numLayers  = [];
numNeurons = cfg.numNeurons;  % <-- ablation controlled

parameters = create_mlp(numLayers, numNeurons, nfeatures, numOutputs);

%% ---------------- Training options ----------------
executionEnvironment = select_execution_environment();

Nepochs              = 1000;
alpha_leaky          = 0.01;

initialLearnRate = 0.01;
decayRate        = 0.0001;

cacheNSteps   = 20;
itersPerEpoch = 200;

averageGrad    = [];
averageSqGrad  = [];
iteration      = 0;

accfun = dlaccelerate(@modelGradients_mlp_MIMOSA_normOnly);

epoch_loss_mean = zeros(Nepochs, 1, 'single');

fprintf('Execution environment: %s\n', char(executionEnvironment));

%% ---------------- Training loop ----------------
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
    cache = cell(cacheNSteps, 2);
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

        dlin  = reshape(dlin3,  [nfeatures,  miniBatchSize * Npt]);
        dlout = reshape(dlout3, [numOutputs, miniBatchSize * Npt]);

        dlin  = dlarray(dlin,  'CB');
        dlout = dlarray(dlout, 'CB');

        [gradients, loss, loss_norm] = dlfeval( ...
            accfun, parameters, dlin, dlout, alpha_leaky); %#ok<NASGU>

        learningRate = initialLearnRate / (1 + decayRate * iteration);

        [parameters, averageGrad, averageSqGrad] = adamupdate( ...
            parameters, gradients, averageGrad, averageSqGrad, iteration, learningRate);

        lossVal = double(gather(extractdata(loss)));
        epoch_loss_total = epoch_loss_total + lossVal;
    end

    epoch_loss_mean(epoch) = single(epoch_loss_total / itersPerEpoch);

    if mod(epoch, 10) == 0 || epoch == 1 || epoch == Nepochs
        fprintf('Run %s | Epoch %d / %d | loss = %.6g | batch = %d\n', ...
            run_name, epoch, Nepochs, epoch_loss_mean(epoch), miniBatchSize);
    end
end

%% ---------------- Final model ----------------
dlnet = struct();
dlnet.parameters = parameters;
dlnet.epoch      = Nepochs;
dlnet.alpha      = alpha_leaky;

%% ---------------- Validation ----------------
Nvalid = 20000;

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

%% ---------------- Save outputs ----------------
save(fullfile(out_dir, [run_name '_final.mat']), ...
    'dlnet', 'epoch_loss_mean', 'valMetrics', 'valPerPoint', 'cfg', 'feature_idx', 'seed');

write_metrics_text(fullfile(out_dir, [run_name '_summary.txt']), ...
    run_name, cfg, feature_idx, epoch_loss_mean(end), valMetrics);

fprintf('\nFinished %s\n', run_name);
fprintf('Val Complex-NRSE = %.6g\n', valMetrics.complex_nrse);
fprintf('Val Complex MAE  = %.6g\n', valMetrics.complex_mae);
fprintf('Val mean|dphi|   = %.6g\n', valMetrics.mean_abs_dphi);
fprintf('Val Mag MAE      = %.6g\n', valMetrics.mag_mae);

end

%% ========================================================================
%% Local functions
%% ========================================================================

function cfg_id = resolve_cfg_id(cfg_name, cfgs)
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

function name_out = sanitize_name(name_in)
% English comments only
name_out = regexprep(char(name_in), '[^a-zA-Z0-9_\-]+', '-');
end

function executionEnvironment = select_execution_environment()
% English comments only
executionEnvironment = "cpu";
try
    nGPU = gpuDeviceCount;
    if nGPU > 0
        gpuDevice([]);
        executionEnvironment = "gpu";
    end
catch
    executionEnvironment = "cpu";
end
end

function write_metrics_text(txtFile, run_name, cfg, feature_idx, final_train_loss, valMetrics)
% English comments only

fid = fopen(txtFile, 'w');
if fid < 0
    warning('Cannot open summary text file for writing: %s', txtFile);
    return;
end

fprintf(fid, 'Run name: %s\n', run_name);
fprintf(fid, 'Config name: %s\n', cfg.name);
fprintf(fid, 'Feature idx: %s\n', mat2str(feature_idx));
fprintf(fid, 'Feature count: %d\n', numel(feature_idx));
fprintf(fid, 'Frequency mode: %s\n', cfg.freq_mode);
fprintf(fid, 'Time mode: %s\n', cfg.time_mode);
fprintf(fid, 'Neuron div: %d\n', cfg.neuron_div);
fprintf(fid, 'numNeurons: %s\n', mat2str(cfg.numNeurons));
fprintf(fid, 'Final train loss: %.10g\n', final_train_loss);

fprintf(fid, '\nValidation metrics\n');
fprintf(fid, 'mean|s| GT      : %.10g\n', valMetrics.mean_abs_gt);
fprintf(fid, 'mean|s| Pred    : %.10g\n', valMetrics.mean_abs_pred);
fprintf(fid, 'mean|delta s|   : %.10g\n', valMetrics.mean_abs_err);
fprintf(fid, 'mean dphi       : %.10g\n', valMetrics.mean_dphi);
fprintf(fid, 'std dphi        : %.10g\n', valMetrics.std_dphi);
fprintf(fid, 'mean |dphi|     : %.10g\n', valMetrics.mean_abs_dphi);
fprintf(fid, 'complex NRSE    : %.10g\n', valMetrics.complex_nrse);
fprintf(fid, 'complex MAE     : %.10g\n', valMetrics.complex_mae);
fprintf(fid, 'NRSE Re         : %.10g\n', valMetrics.nrse_re);
fprintf(fid, 'NRSE Im         : %.10g\n', valMetrics.nrse_im);
fprintf(fid, 'Magnitude NRSE  : %.10g\n', valMetrics.mag_nrse);
fprintf(fid, 'Magnitude MAE   : %.10g\n', valMetrics.mag_mae);

fclose(fid);
end

function cfgs = build_ablation_configs_arch_feat()
% English comments only
% Build 12 configs: 4 neuron divisors x 3 freq modes
% Fixed time set: PT + prepOH
% No FW feature

baseNeurons = [160 240 320 360 480 520 600];
neuronDivs  = [1 2 4 8];

mandatory = 1:9;

freqA = [10 11];
freqB = 12:15;
freqSets  = {freqA, freqB, [freqA freqB]};
freqNames = {'freqRaw', 'phaseTrig', 'freqRaw+phaseTrig'};

timeSet  = [16, 17:20];
timeName = 'PT+prepOH';

cfgs = struct('id', {}, 'name', {}, 'feature_idx', {}, 'nfeatures', {}, ...
              'freq_mode', {}, 'time_mode', {}, 'neuron_div', {}, 'numNeurons', {});

ic = 0;
for idv = 1:numel(neuronDivs)
    div = neuronDivs(idv);

    nn = round(baseNeurons / div);
    nn(nn < 8) = 8;

    for ifr = 1:numel(freqSets)
        ic = ic + 1;

        feat = [mandatory, freqSets{ifr}, timeSet];
        feat = unique(feat, 'stable');

        cfgs(ic).id         = ic; %#ok<AGROW>
        cfgs(ic).freq_mode  = freqNames{ifr}; %#ok<AGROW>
        cfgs(ic).time_mode  = timeName; %#ok<AGROW>
        cfgs(ic).neuron_div = div; %#ok<AGROW>
        cfgs(ic).numNeurons = nn; %#ok<AGROW>
        cfgs(ic).feature_idx = feat; %#ok<AGROW>
        cfgs(ic).nfeatures   = numel(feat); %#ok<AGROW>
        cfgs(ic).name = sprintf('NDIV%d__%s__%s', div, cfgs(ic).freq_mode, cfgs(ic).time_mode); %#ok<AGROW>
    end
end
end

function [gradients, loss, loss_norm] = modelGradients_mlp_MIMOSA_normOnly(parameters, dlXf, dlRf, alpha)
% English comments only

epsA = dlarray(single(1e-4), 'CB');
tau  = dlarray(single(1e-2), 'CB');

if nargin < 4 || isempty(alpha)
    alpha = 0.01;
end

U = mlp_model_leakyRelu(parameters, dlXf, alpha);

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
nfeatures = numel(feature_idx);

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

dlin  = reshape(dlin3,  [nfeatures,  Nvalid * Npt]);
dlout = reshape(dlout3, [numOutputs, Nvalid * Npt]);

dlX = dlarray(dlin,  'CB');
dlR = dlarray(dlout, 'CB');

if executionEnvironment == "gpu"
    dlX = gpuArray(dlX);
    dlR = gpuArray(dlR);
end

U = mlp_model_leakyRelu(parameters, dlX, alpha_leaky);
U = gather(extractdata(U));
R = gather(extractdata(dlR));

pred_reim = reshape(U, [2, Nvalid, Npt]);
gt_reim   = reshape(R, [2, Nvalid, Npt]);

pred_c = squeeze(pred_reim(1,:,:) + 1i * pred_reim(2,:,:));
gt_c   = squeeze(gt_reim(1,:,:)   + 1i * gt_reim(2,:,:));

err_c = pred_c - gt_c;
err_abs = abs(err_c);
gt_abs  = abs(gt_c);
pred_abs = abs(pred_c);

err_abs_vec = err_abs(:);
gt_abs_vec  = gt_abs(:);
pred_abs_vec = pred_abs(:);

tmp = pred_c .* conj(gt_c);
tmp_vec = tmp(:);
dphi = angle(tmp_vec);

valMetrics = struct();

valMetrics.mean_abs_gt   = mean(gt_abs_vec);
valMetrics.mean_abs_pred = mean(pred_abs_vec);
valMetrics.mean_abs_err  = mean(err_abs_vec);

valMetrics.mean_dphi     = mean(dphi);
valMetrics.std_dphi      = std(dphi);
abs_dphi = abs(dphi);
valMetrics.mean_abs_dphi = mean(abs_dphi);

sum_err2 = sum(err_abs_vec .^ 2);
sum_gt2  = sum(gt_abs_vec  .^ 2);
valMetrics.complex_nrse = sqrt(sum_err2) / sqrt(max(sum_gt2, eps('single')));
valMetrics.complex_mae  = mean(err_abs_vec);

pred_re = real(pred_c);
gt_re   = real(gt_c);
pred_im = imag(pred_c);
gt_im   = imag(gt_c);

pred_re_vec = pred_re(:);
gt_re_vec   = gt_re(:);
pred_im_vec = pred_im(:);
gt_im_vec   = gt_im(:);

re_err_vec = pred_re_vec - gt_re_vec;
im_err_vec = pred_im_vec - gt_im_vec;

valMetrics.nrse_re = sqrt(sum(re_err_vec .^ 2)) / sqrt(max(sum(gt_re_vec .^ 2), eps('single')));
valMetrics.nrse_im = sqrt(sum(im_err_vec .^ 2)) / sqrt(max(sum(gt_im_vec .^ 2), eps('single')));

mag_err_vec = pred_abs_vec - gt_abs_vec;
abs_mag_err_vec = abs(mag_err_vec);

valMetrics.mag_nrse = sqrt(sum(mag_err_vec .^ 2)) / sqrt(max(sum(gt_abs_vec .^ 2), eps('single')));
valMetrics.mag_mae  = mean(abs_mag_err_vec);

nrse_pt = zeros(1, Npt, 'single');
for p = 1:Npt
    diff_pt = pred_c(:,p) - gt_c(:,p);
    diff_abs_pt = abs(diff_pt);
    gt_abs_pt = abs(gt_c(:,p));
    e2 = sum(diff_abs_pt .^ 2);
    g2 = sum(gt_abs_pt   .^ 2);
    nrse_pt(p) = sqrt(e2 / max(g2, eps('single')));
end

perPoint = struct();
perPoint.complex_nrse_pt = nrse_pt;
end

function [dlin3, dlout3] = make_minibatch_MIMOSA14_3pool( ...
    feature_idx, batchSize, seq, ...
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

nFull = 22;
featuresFull = zeros(batchSize, Npt, nFull, 'single');
label        = zeros(batchSize, Npt, numOutputs, 'single');

% ---- Sample fractions ----
MWF = single(rand(1, batchSize) * diff(MWF_range) + min(MWF_range));

IEW = zeros(1, batchSize, 'single');
IEW_min = IEW_range(1);
for k = 1:batchSize
    IEW_max_k = min(IEW_range(2), 1 - MWF(k));
    if IEW_max_k < IEW_min
        IEW_max_k = IEW_min;
    end
    IEW(k) = single(rand() * (IEW_max_k - IEW_min) + IEW_min);
end

% ---- Sample T1 and R2 ----
T1_MW_s  = single(rand(1, batchSize) * diff(T1_MW_range_s)  + min(T1_MW_range_s));
T1_IEW_s = single(rand(1, batchSize) * diff(T1_IEW_range_s) + min(T1_IEW_range_s));

R2_MW_Hz  = single(rand(1, batchSize) * diff(R2_MW_range_Hz)  + min(R2_MW_range_Hz));
R2_IEW_Hz = single(rand(1, batchSize) * diff(R2_IEW_range_Hz) + min(R2_IEW_range_Hz));

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
freqMW_ppm  = single(rand(1, batchSize) * diff(freq_MW_range_ppm)  + min(freq_MW_range_ppm));
freqIEW_ppm = single(rand(1, batchSize) * diff(freq_IEW_range_ppm) + min(freq_IEW_range_ppm));
freqBKG_ppm = single(rand(1, batchSize) * diff(freq_BKG_range_ppm) + min(freq_BKG_range_ppm));

b1 = single(rand(1, batchSize) * diff(b1_range) + min(b1_range));

% ---- Pool fractions ----
f_mw  = MWF;
f_iew = IEW;
f_fw  = 1 - MWF - IEW;

% ---- Local frequencies ----
lf_mw_Hz  = freqMW_ppm  * ppm2Hz;
lf_iew_Hz = freqIEW_ppm * ppm2Hz;
lf_bkg_Hz = freqBKG_ppm * ppm2Hz; %#ok<NASGU>

local_mw_Hz  = lf_mw_Hz;
local_iew_Hz = lf_iew_Hz;
local_fw_Hz  = zeros(size(lf_iew_Hz));

% ---- TE list ----
TE_all = [seq.TE_flash seq.TE_flash seq.TE_flash seq.TE_flash seq.TEs(:)'];

% ---- Decay and phase features ----
T2s_MW_decay  = exp(- (R2sMW(:)  * TE_all));
T2s_IEW_decay = exp(- (R2sIEW(:) * TE_all));

phi_mw_tot  = 2*pi * (local_mw_Hz(:)  * TE_all);
phi_iew_tot = 2*pi * (local_iew_Hz(:) * TE_all);

sinPhi_MW  = sin(phi_mw_tot);
cosPhi_MW  = cos(phi_mw_tot);
sinPhi_IEW = sin(phi_iew_tot);
cosPhi_IEW = cos(phi_iew_tot);

% ---- Simulate pool signals ----
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

% ---- Build feature tensors ----
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

t_acq14_s  = single(seq.t_acq14_s(:).');
t_acq_2d   = repmat(t_acq14_s, batchSize, 1);
tmax       = max(t_acq14_s);
t_acq_norm = t_acq_2d ./ (tmax + eps('single')); %#ok<NASGU>

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

% 10..11 raw freq
featuresFull(:,:,10) = Freq_MW_2d;
featuresFull(:,:,11) = Freq_IEW_2d;

% 12..15 trig phase
featuresFull(:,:,12) = single(sinPhi_MW);
featuresFull(:,:,13) = single(cosPhi_MW);
featuresFull(:,:,14) = single(sinPhi_IEW);
featuresFull(:,:,15) = single(cosPhi_IEW);

% 16..20 PT + prepOH
featuresFull(:,:,16) = PT_2d;
prepOH = zeros(batchSize, Npt, 4, 'like', featuresFull);
prepOH(:,1,1) = 1;
prepOH(:,2,2) = 1;
prepOH(:,3,3) = 1;
prepOH(:,4,4) = 1;
featuresFull(:,:,17:20) = prepOH;

% 21 t_acq_norm exists but not used by your new configs
featuresFull(:,:,21) = t_acq_norm;

% 22 FW exists but NOT used by your new configs
featuresFull(:,:,22) = FW_2d;

features = featuresFull(:,:,feature_idx);

dlin3  = permute(features, [3 1 2]);
dlout3 = permute(label,    [3 1 2]);
end