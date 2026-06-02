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
extraData.freqBKG   = pars0.fieldmap / (gpuMIMOSAMWI_A05.gyro*fixed_params.B0); % in ppm
extraData.pini      =  pars0.pini;
extraData.IWF       = pars0.fiew_init;

scale = prctile(PD_slc(:),98);

extraData.PD       = PD_slc/scale;
extraData.MWF       = PV1_slc;
extraData.b1       = b1_slc;

fitting.usingANN = true;   % MLP
fitting.model= 'invivo_1mm';% 'invivo_1mm','invivo_700um'
[out_full]    = objGPU.estimate(y/scale, mask, extraData, fitting);
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
    extraData.freqBKG   = pars0.fieldmap / (gpuMIMOSAMWI_A05.gyro*fixed_params.B0); % in ppm
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
    
    results_3D.freqMW(slc, :, :) = sq(out_ann.final.freqMW * gpuMIMOSAMWI_A05.gyro * fixed_params.B0);
    results_3D.freqIEW(slc, :, :) = sq(out_ann.final.freqIEW * gpuMIMOSAMWI_A05.gyro * fixed_params.B0);
    results_3D.dfreqBKG(slc, :, :) = sq(out_ann.final.dfreqBKG * gpuMIMOSAMWI_A05.gyro * fixed_params.B0);
    results_3D.dpini(slc, :, :) = sq(out_ann.final.dpini);
    
    results_3D.fieldmap_total(slc, :, :) = sq(out_ann.final.dfreqBKG * gpuMIMOSAMWI_A05.gyro * fixed_params.B0 + pars0.fieldmap);
    results_3D.pini_total(slc, :, :) = sq(out_ann.final.dpini + pars0.pini);
    
    results_3D.T1MW(slc, :, :) = sq(out_ann.final.T1MW * 1000);
    results_3D.T1IEW(slc, :, :) = sq(out_ann.final.T1IEW * 1000);
    results_3D.r2MW(slc, :, :) = sq(out_ann.final.r2MW);
    results_3D.r2IEW(slc, :, :) = sq(out_ann.final.r2IEW);
    end
    
end
