% pulseq v1.4.2 is required to be installed first: https://github.com/pulseq/pulseq

sys = mr.opts('MaxGrad', 62, 'GradUnit', 'mT/m', ...
              'MaxSlew', 160, 'SlewUnit', 'T/m/s', ...
              'rfDeadTime', 100e-6, ...
              'rfRingdownTime', 60e-6, ...     
              'adcDeadTime', 40e-6, ...         
              'adcRasterTime', 2e-6, ...        
              'gradRasterTime', 20e-6, ...     
              'blockDurationRaster', 20e-6, ... 
              'B0', 3.0);     % 2.89

sys_lowsar = mr.opts('MaxGrad', 32, 'GradUnit', 'mT/m', ...
              'MaxSlew', 160, 'SlewUnit', 'T/m/s', ...
              'rfDeadTime', 100e-6, ...
              'rfRingdownTime', 60e-6, ...     
              'adcDeadTime', 40e-6, ...         
              'adcRasterTime', 2e-6, ...        
              'gradRasterTime', 20e-6, ...     
              'blockDurationRaster', 20e-6, ... 
              'B0', 3.0);     % 2.89

seq=mr.Sequence(sys);           % Create a new sequence object
%%
% Load ky-kz trajectory data
traj = readmatrix('csv_imports/cplm/fov_224x192_msize_224x192_tl_127_ncal_0_spiral_30_acc_2.5x2.5_nTR215_nph5.txt');

traj = traj(:,[3,2,1]); 
nETL = 127;
nacq = 5;
nTR = size(traj,1)/nETL/nacq;


% Create sequence object and define FOV, resolution, and other high-level parameters
fov=[240 224 192]*1e-3;         % Define FOV and resolution  1 iso
N = [240 224 192];

Nx = N(1);
Ny = N(2);
Nz = N(3);

bandwidth = 340;                % Hz/pixel
flip_angle = 4;                 % degrees
rfSpoilingInc=117;              % RF spoiling increment
dwell = 10e-6;
Tread = dwell*Nx;
os_factor = 1;                  % readout oversampling amount

%--------------------------------------------------------------------------
% Prepare sequence blocks
%--------------------------------------------------------------------------

deltak = 1 ./ fov;
gx = mr.makeTrapezoid('x',sys,'Amplitude',Nx*deltak(1)/Tread,'FlatTime',ceil(Tread/sys.gradRasterTime)*sys.gradRasterTime ,'system',sys);    % readout gradient
gxPre = mr.makeTrapezoid('x',sys,'Area',-gx.area/2,'system',sys_lowsar);        % Gx prewinder
gxSpoil = mr.makeTrapezoid('x',sys_lowsar,'Area',gx.area, 'system',sys_lowsar);         % Gx spoiler

% Make trapezoids for inner loop to save computation

gyPre = mr.makeTrapezoid('y','Area',deltak(2)*(Ny/2),'Duration',mr.calcDuration(gxPre), 'system',sys_lowsar);
gzPre = mr.makeTrapezoid('z','Area',deltak(3)*(Nz/2),'Duration',mr.calcDuration(gxPre), 'system',sys_lowsar);

gyReph = mr.makeTrapezoid('y','Area',deltak(2)*(Ny/2),'Duration',mr.calcDuration(gxSpoil), 'system',sys_lowsar);
gzReph = mr.makeTrapezoid('z','Area',deltak(3)*(Nz/2),'Duration',mr.calcDuration(gxSpoil), 'system',sys_lowsar);

stepsY=((traj(:,2)-1)-Ny/2)/Ny*2;
stepsZ=((traj(:,3)-1)-Nz/2)/Nz*2;

% Create non-selective pulse
rf = mr.makeBlockPulse(flip_angle*pi/180, sys, 'Duration', 0.1e-3);

% Spoilers after t2prep and IR prep
gslSp_t2prep = mr.makeTrapezoid('z','Amplitude',-42.58*8*1e3,'Risetime',0.84e-3,'Duration',0.84e-3+8e-3+0.84e-3,'system',sys); %Amplitude in Hz/m
gslSp_IRprep = mr.makeTrapezoid('z','Amplitude',-42.58*8*1e3,'Risetime',1e-3,'Duration',1e-3+8e-3+1e-3,'system',sys); %Amplitude in Hz/m

% Analog to digital convertersys
adc = mr.makeAdc(Nx * os_factor,'Duration',Tread,'Delay',gx.riseTime,'system',sys);

% Prep pulses
% T2 prep and IR prep pulse are imported from external txt files
% txt file is in mag and phase, while mr.makeArbitraryRf assumes real and imag
t2prep = readmatrix('csv_imports/T2prep.txt');
Re = t2prep(:,1) .* cos(t2prep(:,2));
Im = t2prep(:,1) .* sin(t2prep(:,2));
t2prep_pulse = mr.makeArbitraryRf((Re+Im*1i).', 380.4*pi/180, 'system',sys, 'dwell', 1e-6);

%rf90 = mr.makeBlockPulse(pi/2,sys,'Duration',3e-4);
%rf90_180PhaseOffset = mr.makeBlockPulse(pi/2,sys,'Duration',3e-4,'PhaseOffset',-180*pi/180);

text = readmatrix('csv_imports/rf90.txt');
Re = text(:,1) .* cos(text(:,2));
Im = text(:,1) .* sin(text(:,2));
rf90 = mr.makeArbitraryRf((Re+Im*1i).', pi/2, 'system',sys, 'dwell', 1e-6);
rf90_180PhaseOffset = mr.makeArbitraryRf((Re+Im*1i).', pi/2, 'system',sys, 'dwell', 1e-6, 'PhaseOffset',-180*pi/180);

IRprep = readmatrix('csv_imports/invpulse.txt');
Re = IRprep(:,1) .* cos(IRprep(:,2));
Im = IRprep(:,1) .* sin(IRprep(:,2));
IRprep_pulse = mr.makeArbitraryRf((Re+Im*1i).', 1500*pi/180, 'system',sys, 'dwell', 1e-5);


% #############Multiechoes#############
esp_mte = 2.58;
TEs = [1.98:esp_mte:26.5] * 1e-3;%
TR_mte = 27.7e-3;
nechoes = length(TEs);
gxFlyBack = mr.makeTrapezoid('x','Area',-gx.area,'system',sys);
gxm = mr.scaleGrad(gx,-1);

% gradient spoiling
if mod(length(TEs),2)==0, spSign=-1; else, spSign=1; end
gxSpoil_mgre = mr.makeTrapezoid('x','Area',gx.area*spSign,'system',sys_lowsar);  


delayTE_mte = zeros(nechoes,1);
delayTE_mte(1) = ceil((TEs(1) - (mr.calcDuration(rf) - mr.calcRfCenter(rf)-rf.delay) - mr.calcDuration(gxPre) - mr.calcDuration(gx)/2)/seq.gradRasterTime)*seq.gradRasterTime;
for c = 2:nechoes
    delayTE_mte(c) = ceil(( TEs(c) - TEs(c-1) - mr.calcDuration(gx) )/seq.gradRasterTime)*seq.gradRasterTime;
    if delayTE_mte(c) < 0
        disp(['echo ', num2str(c), ' cannot be fit'])
    else
        disp(['echo ', num2str(c), ' delay ', num2str(1e3*delayTE_mte(c)), ' ms'])
    end
end

delayTR_mte = ceil((TR_mte - mr.calcDuration(rf) - mr.calcDuration(gxPre)  ...
    - mr.calcDuration(gxSpoil) - sum(delayTE_mte) - mr.calcDuration(gx)*length(TEs)...
    )/seq.gradRasterTime)*seq.gradRasterTime;

disp(['delay TR: ', num2str(delayTR_mte*1e3), ' ms'])
dTR=mr.makeDelay(delayTR_mte);
%#################################################

%-------------------------------------------------------------------------
% Adjust sequence timings
%--------------------------------------------------------------------------

esp = 5.8e-3;
gap_between_readouts = 900e-3;

% TE = 20 ms
delay_pre_t2prep_TE20 = 80e-3;
delay_1_t2prep_TE20  =   11e-3 + 80e-6 - 10e-3;      % 11.08 ms
delay_2_t2prep_TE20  =   5e-3;              % 25 ms
delay_3_t2prep_TE20  =   14e-3 - 80e-6 - 10e-3;      % 13.92 ms

% TE = 80 ms
delay_pre_t2prep_TE80 = 20e-3;
delay_1_t2prep_TE80  =   11e-3 + 80e-6 - 2.5e-3;      % 11.08 ms
delay_2_t2prep_TE80  =   20e-3;              % 25 ms
delay_3_t2prep_TE80  =   14e-3 - 80e-6 - 2.5e-3;      % 13.92 ms


delay_IRprep    =   100e-3 - mr.calcDuration(IRprep_pulse)/2;         % gap between end of inversion and start of readout#2 
delay_TE        =   0;
delay_TRinner   =   esp - (mr.calcDuration(rf) + delay_TE + mr.calcDuration(gxPre)+mr.calcDuration(gx)+mr.calcDuration(gxSpoil));         
delay_TRouter = 1e-3;
delT_M3_M4      =   gap_between_readouts - esp*nETL - mr.calcDuration(IRprep_pulse) - delay_IRprep;     % between end of readout#1 and start of inversion
delT_M3_M4      =   delT_M3_M4 - 0.22e-3;
delT_M13_2end   =   53.5e-3;

% %%  calibaration scan
%% %--------------------------------------------------------------------------
Ny_ref = 32;
Nz_ref = 32;

% Calculate timing
TE_ref = 5e-3;
TR_ref = 12e-3;

rf_acs = rf;

delayTE_ref = ceil( (TE_ref - mr.calcDuration(rf_acs) + mr.calcRfCenter(rf_acs) + rf_acs.delay - mr.calcDuration(gxPre)  ...
    - mr.calcDuration(gx)/2)/seq.gradRasterTime)*seq.gradRasterTime;

if delayTE_ref < 0
    disp('error: acs delay TE is negative')    
else
    disp(['acs delay TE: ', num2str(1e3*delayTE_ref), ' ms'])
end


delayTR_ref = ceil((TR_ref - mr.calcDuration(rf_acs) - delayTE_ref- mr.calcDuration(gxPre) ...
    - mr.calcDuration(gx) - mr.calcDuration(gxSpoil))/seq.gradRasterTime)*seq.gradRasterTime;

if delayTR_ref < 0
    disp('error: acs delay TR is negative')    
else
    disp(['acs delay TR: ', num2str(delayTR_ref*1e3), ' ms'])
end


%--------------------------------------------------------------------------
%  dummies
%--------------------------------------------------------------------------
Ndummy_acs = 50;
areaY = ((0:Ny-1)-Ny/2)*deltak(2);
areaZ = ((0:Nz-1)-Nz/2)*deltak(3);
gyPre_acs = mr.makeTrapezoid('y','Area',areaY(floor(Ny/2)),'Duration',mr.calcDuration(gxPre), 'system',sys);
gyReph_acs = mr.makeTrapezoid('y','Area',-areaY(floor(Ny/2)),'Duration',mr.calcDuration(gxPre), 'system',sys);

% Drive magnetization to steady state
rf_phase=0;
rf_inc=0;
for iY = 1:Ndummy_acs
    % RF
    % RF spoiling
    rf_acs.phaseOffset=rf_phase/180*pi;
    % adc.phaseOffset=rf_phase/180*pi;
    rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
    rf_phase=mod(rf_phase+rf_inc, 360.0);       %increment RF phase
    seq.addBlock(rf_acs);


    % Gradients    
    seq.addBlock(gxPre,gyPre_acs);                  % add Gx pre-winder, go to desired ky
    seq.addBlock(mr.makeDelay(delayTE_ref));    % add delay needed before the start of readout

    seq.addBlock(gx);                           % add readout Gx

    seq.addBlock(gyReph_acs,gxSpoil);               % add Gx spoiler, and go back to DC in ky
    seq.addBlock(mr.makeDelay(delayTR_ref))     % add delay to the end of TR
end


temp = 1:Ny;
iY_ref_indices = temp(1+end/2-Ny_ref/2:end/2+Ny_ref/2);

temp = 1:Nz;
iZ_ref_indices = temp(1+end/2-Nz_ref/2:end/2+Nz_ref/2);

% Make trapezoids for inner loop to save computation
for iY = iY_ref_indices
    gyPre_acs(iY) = mr.makeTrapezoid('y','Area',areaY(iY),'Duration',mr.calcDuration(gxPre), 'system',sys);
    gyReph_acs(iY) = mr.makeTrapezoid('y','Area',-areaY(iY),'Duration',mr.calcDuration(gxPre), 'system',sys);
end

%--------------------------------------------------------------------------
% ref data:
%--------------------------------------------------------------------------

mask_ref = zeros([Ny,Nz]);

for iZ = iZ_ref_indices
    % Gz blips to go desired kz, and to come back to DC in kz
    gzPre_acs = mr.makeTrapezoid('z','Area',areaZ(iZ),'Duration',mr.calcDuration(gxPre), 'system',sys);         
    gzReph_acs = mr.makeTrapezoid('z','Area',-areaZ(iZ),'Duration',mr.calcDuration(gxPre), 'system',sys);

    for iY = iY_ref_indices
        % disp(['kZ: ', num2str(iZ), ' kY: ', num2str(iY)])
        mask_ref(iY,iZ) = 1;

        % RF spoiling
        rf_acs.phaseOffset=rf_phase/180*pi;
        adc.phaseOffset=rf_phase/180*pi;
        rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
        rf_phase=mod(rf_phase+rf_inc, 360.0);       %increment RF phase
        % adc.phaseOffset = rf.phaseOffset;

        % Excitation
        seq.addBlock(rf_acs);

        % Encoding
        seq.addBlock(gxPre,gyPre_acs(iY),gzPre_acs);        % Gz, Gy blips, Gx pre-winder
        seq.addBlock(mr.makeDelay(delayTE_ref));    % delay until readout

        seq.addBlock(gx,adc);                       % readout

        seq.addBlock(gyReph_acs(iY),gzReph_acs,gxSpoil);% -Gz, -Gy blips, Gx spoiler
        seq.addBlock(mr.makeDelay(delayTR_ref))     % wait until end of TR
    end
end

mosaic(mask_ref, 1, 1, 20), 



%% ########## Dummy scan
nDummies = 1;
useAdc = 0;% default setting is use ADC
mask_traj = zeros([N(2) N(3) nacq]);% to check mask
for iZ = 1:nDummies
    rf_phase=0;
    rf_inc=0;
    
    % T2 prep pulse TE 20
    seq.addBlock(mr.makeDelay(delay_pre_t2prep_TE20));
    seq.addBlock(rf90,mr.makeDelay(delay_1_t2prep_TE20));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_2_t2prep_TE20));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_2_t2prep_TE20));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_2_t2prep_TE20));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_3_t2prep_TE20));
    seq.addBlock(rf90_180PhaseOffset);
    seq.addBlock(gslSp_t2prep);

    ind_acq = 1;% 
    [rf_phase, rf_inc, mask_traj] = addAcq(seq, nETL, iZ, rf, adc, rfSpoilingInc, rf_phase, rf_inc, stepsZ, stepsY, gxPre, gx, gxSpoil, delay_TE, delay_TRinner, gyPre,gyReph,gzPre,gzReph,nacq,ind_acq,traj,nTR,mask_traj,useAdc);

    % T2 prep pulse TE 80
    seq.addBlock(mr.makeDelay(delay_pre_t2prep_TE80));
    seq.addBlock(rf90,mr.makeDelay(delay_1_t2prep_TE80));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_2_t2prep_TE80));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_2_t2prep_TE80));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_2_t2prep_TE80));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_3_t2prep_TE80));
    seq.addBlock(rf90_180PhaseOffset);
    seq.addBlock(gslSp_t2prep);

    % Contrast 2
    ind_acq = 2;
    [rf_phase, rf_inc, mask_traj] = addAcq(seq, nETL, iZ, rf, adc, rfSpoilingInc, rf_phase, rf_inc, stepsZ, stepsY, gxPre, gx, gxSpoil, delay_TE, delay_TRinner, gyPre,gyReph,gzPre,gzReph,nacq,ind_acq,traj,nTR,mask_traj,useAdc);
    seq.addBlock(mr.makeDelay(delT_M3_M4));

    % IR prep
    seq.addBlock(IRprep_pulse);
    seq.addBlock(gslSp_IRprep,mr.makeDelay(delay_IRprep));

    % Contrast 3
    ind_acq = 3;
    [rf_phase, rf_inc, mask_traj] = addAcq(seq, nETL, iZ, rf, adc, rfSpoilingInc, rf_phase, rf_inc, stepsZ, stepsY, gxPre, gx, gxSpoil, delay_TE, delay_TRinner, gyPre,gyReph,gzPre,gzReph,nacq,ind_acq,traj,nTR,mask_traj,useAdc);
    seq.addBlock(mr.makeDelay(delay_TRouter))

    % Contrast 4
    ind_acq = 4;
    [rf_phase, rf_inc, mask_traj] = addAcq(seq, nETL, iZ, rf, adc, rfSpoilingInc, rf_phase, rf_inc, stepsZ, stepsY, gxPre, gx, gxSpoil, delay_TE, delay_TRinner, gyPre,gyReph,gzPre,gzReph,nacq,ind_acq,traj,nTR,mask_traj,useAdc);
    seq.addBlock(mr.makeDelay(delay_TRouter));

  
    % Contrast 5-14
    ind_acq = 5; 
    [rf_phase, rf_inc, mask_traj] = addAcq_mte(seq, nETL, iZ, rf, adc, rfSpoilingInc, rf_phase, rf_inc, stepsZ, stepsY, gxPre, gx, gxSpoil, delay_TE, delay_TRinner, gyPre,gyReph,gzPre,gzReph,TEs, gxm,gxSpoil_mgre,delayTE_mte, delayTR_mte,nacq,ind_acq,traj,nTR,mask_traj,useAdc);
    seq.addBlock(mr.makeDelay(delT_M13_2end));

end

%% --------------------------------------------------------------------------
% Build sequence
%--------------------------------------------------------------------------
mask_traj = zeros([N(2) N(3) nacq]);% to check mask
useAdc = 1;% default setting is use ADC
for iZ = 1:nTR
    rf_phase=0;
    rf_inc=0;

    % T2 prep pulse TE 20
    seq.addBlock(mr.makeDelay(delay_pre_t2prep_TE20));
    seq.addBlock(rf90,mr.makeDelay(delay_1_t2prep_TE20));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_2_t2prep_TE20));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_2_t2prep_TE20));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_2_t2prep_TE20));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_3_t2prep_TE20));
    seq.addBlock(rf90_180PhaseOffset);
    seq.addBlock(gslSp_t2prep);

    ind_acq = 1;
    [rf_phase, rf_inc, mask_traj] = addAcq(seq, nETL, iZ, rf, adc, rfSpoilingInc, rf_phase, rf_inc, stepsZ, stepsY, gxPre, gx, gxSpoil, delay_TE, delay_TRinner, gyPre,gyReph,gzPre,gzReph,nacq,ind_acq,traj,nTR,mask_traj,useAdc);

    % T2 prep pulse TE 80
    seq.addBlock(mr.makeDelay(delay_pre_t2prep_TE80));
    seq.addBlock(rf90,mr.makeDelay(delay_1_t2prep_TE80));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_2_t2prep_TE80));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_2_t2prep_TE80));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_2_t2prep_TE80));
    seq.addBlock(t2prep_pulse,mr.makeDelay(delay_3_t2prep_TE80));
    seq.addBlock(rf90_180PhaseOffset);
    seq.addBlock(gslSp_t2prep);

    % Contrast 2
    ind_acq = 2;% 
   [rf_phase, rf_inc, mask_traj] = addAcq(seq, nETL, iZ, rf, adc, rfSpoilingInc, rf_phase, rf_inc, stepsZ, stepsY, gxPre, gx, gxSpoil, delay_TE, delay_TRinner, gyPre,gyReph,gzPre,gzReph,nacq,ind_acq,traj,nTR,mask_traj,useAdc);
    seq.addBlock(mr.makeDelay(delT_M3_M4));

    % IR prep
    seq.addBlock(IRprep_pulse);
    seq.addBlock(gslSp_IRprep,mr.makeDelay(delay_IRprep));

    % Contrast 3
    ind_acq = 3;
  [rf_phase, rf_inc, mask_traj] = addAcq(seq, nETL, iZ, rf, adc, rfSpoilingInc, rf_phase, rf_inc, stepsZ, stepsY, gxPre, gx, gxSpoil, delay_TE, delay_TRinner, gyPre,gyReph,gzPre,gzReph,nacq,ind_acq,traj,nTR,mask_traj,useAdc);
    seq.addBlock(mr.makeDelay(delay_TRouter))

    % Contrast 4
    ind_acq = 4;% 
   [rf_phase, rf_inc, mask_traj] = addAcq(seq, nETL, iZ, rf, adc, rfSpoilingInc, rf_phase, rf_inc, stepsZ, stepsY, gxPre, gx, gxSpoil, delay_TE, delay_TRinner, gyPre,gyReph,gzPre,gzReph,nacq,ind_acq,traj,nTR,mask_traj,useAdc);
    seq.addBlock(mr.makeDelay(delay_TRouter));

    % Contrast 5-14
    ind_acq = 5;% 
    [rf_phase, rf_inc, mask_traj] = addAcq_mte(seq, nETL, iZ, rf, adc, rfSpoilingInc, rf_phase, rf_inc, stepsZ, stepsY, gxPre, gx, gxSpoil, delay_TE, delay_TRinner, gyPre,gyReph,gzPre,gzReph,TEs, gxm,gxSpoil_mgre,delayTE_mte, delayTR_mte,nacq,ind_acq,traj,nTR,mask_traj,useAdc);
    seq.addBlock(mr.makeDelay(delT_M13_2end));

end

%--------------------------------------------------------------------------
% Check timing and write sequence
%--------------------------------------------------------------------------

% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

% [pns_ok, pns_n, pns_c, tpns]=seq.calcPNS('MP_GPA_K2309_2250V_951A_AS82.asc'); % prisma
%  
% if (pns_ok)
%      fprintf('PNS check passed successfully\n');
% else
%      fprintf('PNS check failed! The sequence will probably be stopped by the Gradient Watchdog\n');
% end
% 
% gradSpectrum

% Set definitions
seq.setDefinition('FOV', fov);
seq.setDefinition('Matrix', N);
seq.setDefinition('nETL', nETL);
seq.setDefinition('nTR', nTR);
seq.setDefinition('traj_y', traj(:,2));
seq.setDefinition('traj_z', traj(:,3));
seq.setDefinition('os_factor', os_factor);
%########
seq.setDefinition('TES_mte', TEs);
seq.setDefinition('num_echoes', nechoes);
seq.setDefinition('esp_mte', esp_mte);
seq.setDefinition('TR_mte', TR_mte);

% plot
seq.plot('TimeRange',[9.8 10],'timeDisp','ms');

% Write to pulseq file
filename = strrep(mfilename, 'write', '');


%% save

seq.write(['MWF_MIMOSA_1iso','.seq']);

warning('OFF', 'mr:restoreShape')


%%
%--------------------------------------------------------------------------
% Functions
%--------------------------------------------------------------------------



function [rf_phase, rf_inc, mask_traj] = addAcq(seq, nETL, iZ, rf, adc, rfSpoilingInc, rf_phase, rf_inc, stepsZ, stepsY, gxPre, gx, gxSpoil, delay_TE, delay_TRinner, gyPre,gyReph,gzPre,gzReph,nacq,ind_acq,traj,nTR,mask_traj,useAdc)
for iY = 1:nETL


    index = iY+ (ind_acq-1)*nTR*nETL + (iZ-1)*nETL;% $ cplm


     mask_traj(traj(index,2),traj(index,3),ind_acq) = 1;


    % RF spoiling
    rf.phaseOffset=rf_phase/180*pi;
    adc.phaseOffset=rf_phase/180*pi;
    rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
    rf_phase=mod(rf_phase+rf_inc, 360.0);       %increment RF phase

    % Excitation
    seq.addBlock(rf);

    % Encoding
    seq.addBlock(gxPre, ...
        mr.scaleGrad(gyPre,-stepsY(index)), ...
        mr.scaleGrad(gzPre,-stepsZ(index)));    % Gz, Gy blips, Gx pre-winder

%     seq.addBlock(gx, adc);                      % Gx readout
    if useAdc
        seq.addBlock(gx, adc);                      % Gx readout
    else
        seq.addBlock(gx);                           % Gx readout
    end
    
    seq.addBlock(gxSpoil, ...
        mr.scaleGrad(gyReph,stepsY(index)), ...
        mr.scaleGrad(gzReph,stepsZ(index)));    % -Gz, -Gy blips, Gx spoiler

    seq.addBlock(mr.makeDelay(delay_TRinner));  % wait until desired echo spacing
end
end
%##########
function [rf_phase, rf_inc, mask_traj] = addAcq_mte(seq, nETL, iZ, rf, adc, rfSpoilingInc, rf_phase, rf_inc, stepsZ, stepsY, gxPre, gx, gxSpoil, delay_TE, delay_TRinner, gyPre,gyReph,gzPre,gzReph,TEs, gxm,gxSpoil_mgre,delayTE_mte, delayTR_mte,nacq,ind_acq,traj,nTR,mask_traj,useAdc)
for iY = 1:nETL

    % Calculate index for ky-kz look up table

    % RF spoiling
    rf.phaseOffset=rf_phase/180*pi;
    adc.phaseOffset=rf_phase/180*pi;
    rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
    rf_phase=mod(rf_phase+rf_inc, 360.0);       %increment RF phase

    % Excitation
    seq.addBlock(rf);
    
    % multiecho Encoding
    for c=1:length(TEs) % loop over TEs

        %###########
        % Calculate index for ky-kz look up table
        ind_acq_me = ind_acq;

        index = (ind_acq_me-1)*nTR*nETL - iY + 1 + iZ*nETL;
    
        % disp and show
%         disp(['acq_' num2str(ind_acq_me) ',index=',num2str(index)]);
%         figure(4);hold on;scatter(traj(index,2),traj(index,3),10,ind_acq_me);
        mask_traj(traj(index,2),traj(index,3),ind_acq_me) = 1;

        if c==1 
            seq.addBlock(gxPre, ...
                mr.scaleGrad(gyPre,-stepsY(index)), ...
                mr.scaleGrad(gzPre,-stepsZ(index)));    % Gz, Gy blips, Gx pre-winder
            seq.addBlock(mr.makeDelay(delayTE_mte(c)));
        else
            seq.addBlock(mr.makeDelay(delayTE_mte(c)));

        end
        if useAdc
            if mod(c,2)==0
                seq.addBlock(gxm, adc); 
            else
                seq.addBlock(gx, adc); 
            end
                                 % Gx readout
        else
            if mod(c,2)==0
                seq.addBlock(gxm); 
            else
                seq.addBlock(gx); 
            end                          % Gx readout
        end
        
        ind_pre = index;
    end

    seq.addBlock(gxSpoil_mgre, ...
        mr.scaleGrad(gyReph,stepsY(index)), ...
        mr.scaleGrad(gzReph,stepsZ(index)));    % -Gz, -Gy blips, Gx spoiler

    seq.addBlock(mr.makeDelay(delayTR_mte));  % wait until desired echo spacing
end
end