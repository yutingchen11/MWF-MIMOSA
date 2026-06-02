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
% qsm = Data.qsm_map;
x_tot = Data.x_tot;

