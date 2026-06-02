classdef gpuMIMOSAMWI_A05 < handle

    properties (Constant)
            gyro = 42.57747892;
    end

    % here define the fitting parameters and their conditions
    properties
        % default model parameters and estimation boundary
        % S0        : T1w signal [a.u.] 
        % MWF       : myelin water fraction [0,1]
        % IEW       : intracellular volume ratio (=Vic or ICVF in DWI) [0,1]
        % r2sMW     : R2* MW [1/s]
        % r2sIEW     : R2* IW [1/s] 
        % r2sIEW     : R2* EW [1/s] 
        % freqMW    : frequency MW [ppm], range: -6:30 Hz
        % freqIEW    : frequency IW [ppm], range: -12:6 Hz
        % dfreqBKG  : background frequency in addition to the one provided [ppm]
        % dpini     : B1 phase offset in addition to the one provided [rad]
        modelParams     = {'S0' ;   'MWF';  'IEW';  'r2sMW';'r2sIEW';  'freqMW'; 'freqIEW'; 'dfreqBKG'; 'dpini';  'T1MW'; 'T1IEW';  'r2MW';'r2IEW'}; % modelParams;
        ub              = [    1;     0.3;      1;      200;      50;      0.25;      0.05;        0.4;    pi/2;   0.4;     1.5;        100;    17; ];
        lb              = [ 1e-8;    1e-8;   1e-8;       50;       2;     -0.05;      -0.1;       -0.4;   -pi/2;  0.1;     0.7;        25;    10; ];
        % ub              = [    10;     0.3;      1;      200;      50;      0.25*1.5;      0.05*1.5;        0.4*1.5;    pi];% yc: increase bak range
        % lb              = [ 1e-8;    1e-8;   1e-8;       50;       2;     -0.05*1.5;      -0.1*1.5;        -0.4*1.5;   -pi];
       
        % % startPoint      = [    0.1;     0.1;    0.6;      100;      15;      0.04;         0;          0;       0];
        % ub              = [    10;     0.3;      1;      200;      66;      0.25;      0.05;        0.4;    pi];% yc: -pi 2 pi
        % lb              = [ 1e-8;    1e-8;   1e-8;       66;       14;     -0.05;      -0.1;       -0.4;   -pi];
        startPoint      = [    0.1;     0.1;    0.8;      100;      15;      0.04;         0;          0;       0;  0.15;    1;         67;     14];

    end

    properties (GetAccess = public, SetAccess = protected)

        % B0      = 3;            % T
        B0      = 2.89;            % T, for siemens, yc
        x_i     = -0.1;         % ppm
        x_a     = -0.1;         % ppm
        E       = 0.02;         % ppm, exchange induced frequency shift
        rho_mw  = 0.36/0.86;    % ratio, myelin water proton ratio
        B0dir   = [0;0;1];      % unit vector [x,y,z]

        thres_R2star = 2;

        te
        
    end
    
    methods

        % constructuor
        function this = gpuMIMOSAMWI_A05(te,fixed_params)
        % MWF_MIMOSA 3-pool model
        % Output
        % ----------
        % this      : object of a fitting class

            this.te     = single(te(:));
            % fixed tissue and scanner parameters
            if nargin == 2
                if isfield(fixed_params,'x_i');         this.x_i            = single(fixed_params.x_i);             end
                if isfield(fixed_params,'x_a');         this.x_a            = single(fixed_params.x_a);             end
                if isfield(fixed_params,'E');           this.E              = single(fixed_params.E);               end
                if isfield(fixed_params,'rho_mw');      this.rho_mw         = single(fixed_params.rho_mw);          end
                if isfield(fixed_params,'B0');          this.B0             = single(fixed_params.B0);              end
                if isfield(fixed_params,'B0dir');       this.B0dir          = single(fixed_params.B0dir);           end
                if isfield(fixed_params,'thres_R2s');   this.thres_R2star   = single(fixed_params.thres_R2star);    end
            end
        end
        
        % update properties according to lmax
        function this = updateProperty(this, fitting)

            if fitting.isComplex == 0
                for kpar = {'dfreqBKG','dpini'}
                    idx = find(ismember(this.modelParams,kpar));
                    this.modelParams(idx)    = [];
                    this.lb(idx)              = [];
                    this.ub(idx)              = [];
                    this.startPoint(idx)      = [];
                end
            end

            % DIMWI
            if fitting.DIMWI.isFitFreqIEW == 0
                idx = find(ismember(this.modelParams,'freqIEW'));
                this.modelParams(idx)    = [];
                this.lb(idx)              = [];
                this.ub(idx)              = [];
                this.startPoint(idx)      = [];
            end

            if fitting.DIMWI.isFitFreqMW == 0
                idx = find(ismember(this.modelParams,'freqMW'));
                this.modelParams(idx)    = [];
                this.lb(idx)              = [];
                this.ub(idx)              = [];
                this.startPoint(idx)      = [];
            end

            if fitting.DIMWI.isFitIEW == 0
                idx = find(ismember(this.modelParams,'IEW'));
                this.modelParams(idx)    = [];
                this.lb(idx)              = [];
                this.ub(idx)              = [];
                this.startPoint(idx)      = [];
            end

            if fitting.DIMWI.isFitR2sIEW == 0
                idx = find(ismember(this.modelParams,'r2sIEW'));
                this.modelParams(idx)    = [];
                this.lb(idx)              = [];
                this.ub(idx)              = [];
                this.startPoint(idx)      = [];
            end

        end

        % display some info about the input data and model parameters
        function display_data_model_info(this)

            disp('================================');
            disp('MWF-MIMOSA with askAdam solver');
            disp('================================');
            
            disp('----------------')
            disp('Data Information');
            disp('----------------')
            disp([  'Field strength (T)                     : ' num2str(this.B0)]);
            fprintf('Echo time, TE (ms)                     : [%s] \n',num2str((this.te*1e3).',' %.2f'));
            
            disp('---------------------')
            disp('Parameter to be fixed')
            disp('---------------------')
            disp(['Relative myelin water density            : ' num2str(this.rho_mw)]);
            disp(['Myelin isotropic susceptibility (ppm)    : ' num2str(this.x_i)]);
            disp(['Myelin anisotropic susceptibility (ppm)  : ' num2str(this.x_a)]);
            disp(['Exchange term (ppm)                      : ' num2str(this.E)]);
            disp('---------------------')

        end

        %% higher-level data fitting functions
        % Wrapper function of fit to handle image data; automatically segment data and fitting in case the data cannot fit in the GPU in one go
        function  [out] = estimate(this, data, mask, extraData, fitting)
        % Perform MWF_MIMOSA model parameter estimation based on askAdam
        % Input data are expected in multi-dimensional image
        % 
        % Input
        % -----------
        % data      : 4D multi-echo GRE, [x,y,z,te]
        % mask      : 3D signal mask, [x,y,z]
        % extradata : Optional additional data
        %   .freqBKG: 3D initial estimation of total field [Hz] (highly recommended)
        %   .pini   : 3D initial estimation of B1 offset [rad]  (highly recommended)
        %   .ff     : 3D/4D fibre fraction map, [x,y,z,nF] (for GRE-DIMWI only)
        %   .theta  : 3D/4D angle between B0 and fibre orientation, [x,y,z, nF] (for GRE-DIMWI only)
        %   .IEW    : 3D volume fractino IC/(IC+EC), [x,y,z] (for GRE-DIMWI only)
        % fitting   : fitting algorithm parameters (see fit function)
        % 
        % Output
        % -----------
        % out       : output structure contains all estimation results
        % 
            
           % display basic info
            this.display_data_model_info;

            % get all fitting algorithm parameters 
            fitting = this.check_set_default(fitting,data);

            % get matrix size
            dims = size(data,1:3);

            % make sure input data are valid
            [extraData,mask] = this.validate_data(data,extraData,mask,fitting);

            % compute rotationally invariant signal if needed
            [data, scaleFactor] = this.prepare_data(data,extraData,mask);% yc: nomalized based on input PD

            % mask sure no nan or inf
            [data,mask] = utils.remove_img_naninf(data,mask);

            % convert datatype to single
            data    = single(data);
            mask    = mask >0;

            % determine if we need to divide the data to fit in GPU
            g = gpuDevice; reset(g);
            memoryFixPerVoxel       = 0.0001;   % get this number based on mdl fit
            memoryDynamicPerVoxel   = 0.0001;     % get this number based on mdl fit
            [NSegment,maxSlice]     = utils.find_optimal_divide(mask,memoryFixPerVoxel,memoryDynamicPerVoxel);

            % parameter estimation
            out = [];
            for ks = 1:NSegment

                fprintf('Running #Segment = %d/%d \n',ks,NSegment);
                disp   ('------------------------')
    
                % determine slice# given a segment
                if ks ~= NSegment
                    slice = 1+(ks-1)*maxSlice : ks*maxSlice;
                else
                    slice = 1+(ks-1)*maxSlice : dims(3);
                end
                
                % divide the data
                dwi_tmp     = data(:,:,slice,:);
                mask_tmp    = mask(:,:,slice);
                fields      = fieldnames(extraData); 
                for kfield = 1:numel(fields); extraData_tmp.(fields{kfield}) = single(extraData.(fields{kfield})(:,:,slice,:,:)); end

                % run fitting
                [out_tmp]    = this.fit(dwi_tmp,mask_tmp,fitting,extraData_tmp);

                % restore 'out' structure from segment
                out = utils.restore_segment_structure(out,out_tmp,slice,ks);

            end
            out.mask = mask;
            % rescale S0
            out.final.S0    = out.final.S0 *scaleFactor;
            out.min.S0      = out.min.S0 *scaleFactor;

            % save the estimation results if the output filename is provided
            askadam.save_askadam_output(fitting.outputFilename,out)

        end

        % Data fitting function, can be 2D (voxel) or 4D (image-based)
        function [out] = fit(this,data,mask,fitting,extraData)
        %
        % Input
        % -----------
        % dwi       : S0 normalised 4D dwi images, [x,y,slice,diffusion], 4th dimension corresponding to [Sl0_b1,Sl0_b2,Sl2_b1,Sl2_b2, etc.]; the order of bval must match the order in the constructor gpuNEXI
        % mask      : 3D signal mask, [x,y,slice]
        % fitting   : fitting algorithm parameters
        %   .Nepoch             : no. of maximum iterations, default = 4000
        %   .initialLearnRate   : initial gradient step size, defaulr = 0.01
        %   .decayRate          : decay rate of gradient step size; learningRate = initialLearnRate / (1+decayRate*epoch), default = 0.0005
        %   .convergenceValue   : convergence tolerance, based on the slope of last 'convergenceWindow' data points on loss, default = 1e-8
        %   .convergenceWindow  : number of data points to check convergence, default = 20
        %   .tol                : stop criteria on metric value, default = 1e-3
        %   .lambda             : regularisation parameter, default = 0 (no regularisation)
        %   .TVmode             : mode for TV regulariation, '2D'|'3D', default = '2D'
        %   .regmap             : parameter map used for regularisation, 'fa'|'ra'|'Da'|'De', default = 'fa'
        %   .lossFunction       : loss for data fidelity term, 'L1'|'L2'|'MSE', default = 'L1'
        %   .display            : online display the fitting process on figure, true|false, defualt = false
        %   .isWeighted         : is cost function weighted, true|false, default = true
        %   .weightMethod       : Weighting method, '1stecho'|'norm', default = '1stecho'
        %   .weightPower        : power order of the weight, default = 2
        %   .DIMWI.isFitIEW     : Vic is a free parameter, default = true
        %   .DIMWI.isFitFreqMW  : MW frequency is a free parameter, default = true
        %   .DIMWI.isFitFreqIEW  : IW frequency is a free parameter, default = true
        %   .DIMWI.isFitR2sIEW   : EW R2* is a free parameter, default = true
        % 
        % Output
        % -----------
        % out       : output structure
        %   .final      : final results (see properties for other parameters)
        %       .loss       : final loss metric
        %   .min        : results with the minimum loss metric across all iterations
        %       .loss       : loss metric      
        %
        % Description: askAdam Image-based NEXI model fitting
        %
        % Kwok-Shing Chan @ MGH
        % kchan2@mgh.harvard.edu
        % Date created: 19 July 2024
        % Date modified:
        %
        %
            
            % check GPU
            gpool = gpuDevice;
            
            % get image size
            dims = size(data,1:3);

            %%%%%%%%%%%%%%%%%%%% 1. Validate and parse input %%%%%%%%%%%%%%%%%%%%
            if nargin < 3 || isempty(mask); mask = ones(dims,'logical'); end % if no mask input then fit everthing
            if nargin < 4; fitting = struct(); end

            % get all fitting algorithm parameters 
            fitting             = this.check_set_default(fitting,data);
            % determine fitting parameters
            this                = this.updateProperty(fitting);
            fitting.modelParams = this.modelParams;
            % set fitting boundary if no input from user
            if isempty( fitting.ub); fitting.ub = this.ub(1:numel(this.modelParams)); end
            if isempty( fitting.lb); fitting.lb = this.lb(1:numel(this.modelParams)); end
            
            % set initial starting points
            pars0 = this.determine_x0(data,extraData,fitting);
            % pars0 = this.estimate_prior(data);
            
            %%%%%%%%%%%%%%%%%%%% End 1 %%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%% 2. Setting up all necessary data, run askadam and get all output %%%%%%%%%%%%%%%%%%%%
            % 2.1 setup fitting weights
            w = this.compute_optimisation_weights(data,fitting); % This is a customised funtion

            % split data into real and imaginary parts for complex-valued data
            if fitting.isComplex; data = cat(5,real(data),imag(data)); end

            % 2.2 display optimisation algorithm parameters
            this.display_algorithm_info(fitting)

            % 3. askAdam optimisation main
            askadamObj  = askadam();
            % % mask out data to reduce memory load
            % data = utils.vectorise_NDto2D(data,mask).';
            % if ~isempty(w); w = utils.vectorise_NDto2D(w,mask).'; end
            % fieldname = fieldnames(extraData); for km = 1:numel(fieldname); extraData.(fieldname{km}) = gpuArray(single( utils.vectorise_NDto2D(extraData.(fieldname{km}),mask) ).'); end
            extraData   = utils.masking_ND2GD_preserve_struct(extraData,mask) ;
            out         = askadamObj.optimisation(data, mask, w, pars0, fitting, @this.FWD, fitting, extraData);

            iewRel_final = out.final.IEW;
            mwf_final    = out.final.MWF;
            
            % out.final.IEWrel = iewRel_final;
            out.final.IEW    = (1 - mwf_final) .* iewRel_final;
            out.final.FW     = (1 - mwf_final) .* (1 - iewRel_final);
            
            iewRel_min = out.min.IEW;
            mwf_min    = out.min.MWF;
            
            % out.min.IEWrel = iewRel_min;
            out.min.IEW    = (1 - mwf_min) .* iewRel_min;
            out.min.FW     = (1 - mwf_min) .* (1 - iewRel_min);

            %%%%%%%%%%%%%%%%%%%% End 2 %%%%%%%%%%%%%%%%%%%%

            disp('The process is completed.')
            
            % clear GPU
            reset(gpool)
            
        end

        %% Prior estimation related functions

        % determine how the starting points will be set up
        function x0 = determine_x0(this,y,extraData,fitting) 

            disp('---------------');
            disp('Starting points');
            disp('---------------');

            dims = size(y,1:3);

            if ischar(fitting.start)
                switch lower(fitting.start)
                    case 'prior'
                        % using maximum likelihood method to estimate starting points
                        x0 = this.estimate_prior(y);
    
                    case 'default'
                        % use fixed points
                        fprintf('Using default starting points for all voxels at [%s]: [%s]\n', cell2str(this.modelParams),replace(num2str(this.startPoint(:).',' %.2f'),' ',','));
                        x0 = utils.initialise_x0(dims,this.modelParams,this.startPoint);
                    case 'prior_mimosa' % yc , add
                        x0 = utils.initialise_x0(dims,this.modelParams,this.startPoint);
                    
                        x0.MWF = extraData.MWF;
                        x0.IEW = extraData.IWF;
                        x0.S0 = extraData.PD;
                        % global constant
                        % x0.T1MW = 234e-3;
                        if isfield(extraData, 'freqMW')
                            x0.freqMW = extraData.freqMW;
                        end

                end
            else
                % user defined starting point
                x0 = fitting.start(:);
                fprintf('Using user-defined starting points for all voxels at [%s]: [%s]\n',cell2str(this.modelParams),replace(num2str(x0(:).',' %.2f'),' ',','));
                x0 = utils.initialise_x0(dims,this.modelParams,this.startPoint);

            end
            
            % make sure the input is bounded
            x0 = askadam.set_boundary(x0,fitting.ub,fitting.lb);

            fprintf('Estimation lower bound [%s]: [%s]\n',      cell2str(this.modelParams),replace(num2str(fitting.lb(:).',' %.2f'),' ',','));
            fprintf('Estimation upper bound [%s]: [%s]\n',      cell2str(this.modelParams),replace(num2str(fitting.ub(:).',' %.2f'),'  ',','));
            ('---------------');
        end

        % using maximum likelihood method to estimate starting points
        function pars0 = estimate_prior(this,data)
        % Estimation starting points 

            % data = zeros(size(data));
            % for kt = 1:length(this.te)
            %     data(:,:,:,kt) = smooth3(data(:,:,:,kt));
            % end

            dims = size(data,1:3);

            % initiate starting point of all parameters
            pars0 = utils.initialise_x0(dims,this.modelParams,this.startPoint);
            
            disp('Estimate starting points based on hybrid fixed points/prior information ...')

            % S0
            % [R2s,S0]  =
            % this.R2star_trapezoidal(abs(data(5:end,:)),this.te);% yc, for GD format
            [R2s,S0]  = this.R2star_trapezoidal(abs(data(:,:,:,5:end)),this.te);% yc, for ND format
            S0(isnan(S0)) = 0; S0(isinf(S0)) = 0; S0(S0<0) = 0;
            pars0.S0 = single(S0);

            % R2*IEW
            idx = find(ismember(this.modelParams,'r2sIEW'));
            r2sIEW = R2s - 3;% yc: not sure
            r2sIEW(isnan(r2sIEW)) = single(this.startPoint(idx)); r2sIEW(isinf(r2sIEW)) = single(this.startPoint(idx)); 
            r2sIEW(r2sIEW < this.lb(idx)) = single(this.lb(idx)); r2sIEW(r2sIEW > this.ub(idx)) = single(this.ub(idx));
            pars0.r2sIEW = single(r2sIEW);

            % % R2*EW
            % idx = find(ismember(this.modelParams,'r2sIEW'));
            % if ~isempty(idx)
            %     % if R2*EW is a free parameter then set it
            %     r2sIEW = R2s + 3;
            %     r2sIEW(isnan(r2sIEW)) = single(this.startPoint(idx)); r2sIEW(isinf(r2sIEW)) = single(this.startPoint(idx)); 
            %     r2sIEW(r2sIEW < this.lb(idx)) = single(this.lb(idx)); r2sIEW(r2sIEW > this.ub(idx)) = single(this.ub(idx));
            %     pars0.r2sIEW = single(r2sIEW);
            % % else
            % %     % if R2*EW is not a free parameter then reset R2*EW
            % %     idx = find(ismember(this.modelParams,'r2sIEW'));
            % %     pars0.(this.modelParams{idx}) = single(this.startpoint(idx)*ones(dims));
            % end

            % [~,mwf] = this.superfast_mwi_2m_standard(abs(data),this.te,[]);
            % mwf(mwf>0.15)                   = 0.15;         
            % mwf(and(mwf>=0.05,mwf<=0.1))    = 0.1;   
            % mwf(mwf<0.015)                  = 0.03;
            % pars0.(this.modelParams{2})    = single(mwf);
            

        end

        %% Signal related functions

        % Forward model to generate MWF_MIMOSA signal
        function [s] = FWD(this, pars, fitting, extraData)

            % TE = gpuArray(dlarray( permute(this.te, [2 3 4 1] )));              % TE always on 4th dim

            S0   = pars.S0;
            mwf  = pars.MWF;
            if fitting.DIMWI.isFitIEW; iew  = pars.IEW; else; iew = extraData.IEW; end
            r2sMW   = pars.r2sMW;
            r2sIEW   = pars.r2sIEW;

            if fitting.DIMWI.isFitR2sIEW;    r2sIEW   = pars.r2sIEW;   end
            if fitting.DIMWI.isFitFreqMW;   freqMW  = pars.freqMW;  end
            if fitting.DIMWI.isFitFreqIEW;   freqIEW  = pars.freqIEW;  end
            % external effects
            if ~fitting.isComplex % magnitude fitting
                freqBKG = 0;                          
                pini    = 0;
            else    % other fittings
                freqBKG = pars.dfreqBKG + extraData.freqBKG; 
                pini    = pars.dpini + extraData.pini;
            end
        
            %%%%%%%%%%%%%%%%%%%% Compartmental Signals %%%%%%%%%%%%%%%%%%%%
            % S0MW = S0 .* mwf;
            % S0IW = S0 .* (1-mwf) .* iew;
            % S0CSF = S0 .* (1-mwf) .* (1-iew);
            
            % yc
            % fw = max(1-mwf-iew,0);
            % S0MW = S0 .* mwf;
            % S0IW = S0 .* iew;
            % S0CSF = S0 .* fw;

            % % English comment: enforce 3-compartment simplex (nonnegative + sum-to-one)
            % mwf = max(mwf, 0);
            % mwf = min(mwf, 0.3);
            % 
            % fw  = 1 - mwf - iew;       % guaranteed >= 0
            % fw  = max(fw, 0);          % safety against rounding
            % 
            % iew = max(iew, 0);
            % iew = min(iew, 1 - mwf - fw);   % critical coupling constraint
            % 
            % 
            % 
            % % English comment: use NONNEGATIVE fractions for amplitudes
            % S0MW  = S0 .* mwf;
            % S0IW  = S0 .* iew;
            % S0CSF = S0 .* fw;
            mwf = pars.MWF;
            q_iew = pars.IEW;  % fitted parameter: IEW fraction within non-myelin water
            
            fw  = (1 - mwf) .* (1 - q_iew);
            iew = (1 - mwf) .* q_iew;
            
            S0MW  = S0 .* mwf;
            S0IW  = S0 .* iew;
            S0CSF = S0 .* fw;

            

            %%%%%%%%%%%%%%%%%%%% DIMWI related operations %%%%%%%%%%%%%%%%%%%%
            
            % if use HCFM to derive either freqMW|freqIEW|R2*EW, then computee g-ratio
            if ~fitting.DIMWI.isFitFreqMW || ~fitting.DIMWI.isFitFreqIEW || ~fitting.DIMWI.isFitR2sIEW
                hcfm_obj = HCFM(this.te,this.B0);

                % g-ratio 
                g = hcfm_obj.gratio(abs(S0IW),abs(S0MW)/this.rho_mw);

            end
            
            % extra decay on extracellular water estimated by HCFM 
            if ~fitting.DIMWI.isFitR2sIEW
                
                % assume extracellular water has the same T2* as intra-axonal water
                r2sIEW   = r2sIEW;
                % fibre volume fraction
                fvf     = hcfm_obj.FibreVolumeFraction(abs(S0IW),abs(S0CSF),abs(S0MW)/this.rho_mw);

                % signal dephase in extracellular water due to myelin sheath, Eq.[A7]
                decayEW = hcfm_obj.DephasingExtraaxonal(fvf,g,this.x_i,this.x_a,extraData.theta);

            else
                decayEW = 0;
            end

            % compute frequency shifts given theta
            if ~fitting.DIMWI.isFitFreqMW 

                % in ppm
                freqMW = hcfm_obj.FrequencyMyelin(this.x_i,this.x_a,g,extraData.theta,this.E) / (this.B0*this.gyro);

            end
            if ~fitting.DIMWI.isFitFreqIEW 

                % in ppm
                freqIEW = hcfm_obj.FrequencyAxon(this.x_a,g,extraData.theta) / (this.B0*this.gyro);

            end

            %MIMSOA signal model
            % t2mw = 15e-3;% 25e-3    
            % t2iew =70e-3; t2csf = 500e-3;
            % t1mw = 150e-3; t1iew = 1000e-3; t1csf = 4500e-3;
            t2mw =  1./(pars.r2MW);t2iew =1./(pars.r2IEW); t2csf = 500e-3;
            t1mw = pars.T1MW; t1iew = pars.T1IEW; t1csf = 4500e-3;
            % t1mw = 150e-3; t1iew = 1500e-3; t1csf = 4500e-3;
            % t2smw = 1./(r2sMW + 1/t2mw); t2siew = 1./(r2sIEW + 1/t2iew); t2scsf = 500e-3;
             t2smw = 1./(r2sMW); t2siew = 1./(r2sIEW);t2scsf = 500e-3;
            %%%%%%%%%% simulate signal based on parameter input %%%%%%%%%%%%%%%%%%%%%%%%
            param = fitting.param;

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
            
            prototype = pars.S0; 
            num_vox = numel(pars.S0);

            freq_fw = zeros(num_vox, 1, 'like', prototype);% FW freq

            if size(t1mw, 2) == 1
    
                ones_dl = ones(num_vox, 1, 'like', prototype);
                
                t1mw_dl = ones_dl .* t1mw(:);
                t2mw_dl = ones_dl .* t2mw(:);

                t2smw_dl = reshape(t2smw, [], 1);
                S0_dl    = reshape(S0, [], 1);
                

                tissue_mw = cat(2, t1mw_dl, t2mw_dl, t2smw_dl, S0_dl);
                

                t1iew_dl = ones_dl .* t1iew(:);
                t2iew_dl = ones_dl .* t2iew(:);
                t2siew_dl = reshape(t2siew, [], 1);
                tissue_aew = cat(2, t1iew_dl, t2iew_dl, t2siew_dl, S0_dl);
            
                t1csf_dl = ones_dl .* t1csf;
                t2csf_dl = ones_dl .* t2csf;
                t2scsf_dl = ones_dl .* t2scsf; 
                tissue_fw = cat(2, t1csf_dl, t2csf_dl, t2scsf_dl, S0_dl);
                else
                tissue_mw = cat(2, dlarray(t1mw(:)), dlarray(t2mw(:)), t2smw(:));
                tissue_aew = cat(2, dlarray(t1iew(:)), dlarray(t2iew(:)), t2siew(:));
                tissue_fw = cat(2, dlarray(t1csf(:)), dlarray(t2csf(:)), t2scsf(:));
                
            end



            b1_val = extraData.b1;


            if ~fitting.usingANN
                % inv_eff = 0.8;
                % gen siganl of mw
                % [Mz, Mxy] = sim_mwf_MIMOSA_bp(TR, alpha_deg, esp, turbo_factor, tissue_mw(:,1), tissue_mw(:,2), num_reps, echo2use, TR_mte,esp_mte,TEs,tissue_mw(:,3), gap_between_readouts, time2relax_at_the_end,b1_val, inv_eff,tissue_mw(:,4));
                [Mz, Mxy, Mxy0] = sim_mwf_MIMOSA_bp_mxy0_cx_v2_dlarray(TR, alpha_deg, esp, turbo_factor, tissue_mw(:,1), tissue_mw(:,2), num_reps, echo2use, TR_mte,esp_mte,TEs,tissue_mw(:,3), gap_between_readouts, time2relax_at_the_end,b1_val(:),0.5022,(freqMW(:)).*this.B0.*this.gyro,S0MW(:));              
                % temp = abs(Mxy(:,:,end).');
                % temp1 = abs(Mxy0(:,:,end).');
                % for n = 1:size(temp,1)
                %    temp(n,:)      = temp(n,:) / sum(abs(temp(n,:)).^2)^0.5;
                %    temp1(n,:)      = temp1(n,:) / sum(abs(temp1(n,:)).^2)^0.5;
                % end
      
                temp = Mxy(:,:,end).';
                temp1 = Mxy0(:,:,end).';
                img_sim_mw = temp;
                Smw_mgre = temp1;

                
                % img_sim_mw = reshape(signal,[N(1),N(2),1,14]);
                
                % mosaic( squeeze( img_sim_mw), 2, 5, 1, '',[0 1],0 );colormap gray
                
                
                % gen signal of aew
                [Mz, Mxy, Mxy0] = sim_mwf_MIMOSA_bp_mxy0_cx_v2_dlarray(TR, alpha_deg, esp, turbo_factor, tissue_aew(:,1), tissue_aew(:,2), num_reps, echo2use, TR_mte,esp_mte,TEs,tissue_aew(:,3), gap_between_readouts, time2relax_at_the_end,b1_val(:), 0.8567,(freqIEW(:)).*this.B0.*this.gyro,S0IW(:));
                          
                % temp = abs(Mxy(:,:,end).');
                % temp1 = abs(Mxy0(:,:,end).');
                % for n = 1:size(temp,1)
                %    temp(n,:)      = temp(n,:) / sum(abs(temp(n,:)).^2)^0.5;
                %    temp1(n,:)      = temp1(n,:) / sum(abs(temp1(n,:)).^2)^0.5;
                % end
                temp = Mxy(:,:,end).';
                temp1 = Mxy0(:,:,end).';
                img_sim_aew = temp;
                Siew_mgre = temp1;
                % img_sim_aew = reshape(signal,[N(1),N(2),1,14]);
                % mosaic( squeeze( img_sim_aew), 2, 5, 2, '',[0 1],0 );colormap gray
                
                
                % gen signal of fw
                [Mz, Mxy, Mxy0] = sim_mwf_MIMOSA_bp_mxy0_cx_v2_dlarray(TR, alpha_deg, esp, turbo_factor, tissue_fw(:,1), tissue_fw(:,2), num_reps, echo2use, TR_mte,esp_mte,TEs,tissue_fw(:,3), gap_between_readouts, time2relax_at_the_end,b1_val(:),0.9638,(freq_fw(:)).*this.B0.*this.gyro,S0CSF(:));           
                
                % temp = abs(Mxy(:,:,end).');
                % temp1 = abs(Mxy0(:,:,end).');
                % for n = 1:size(temp,1)
                %    temp(n,:)      = temp(n,:) / sum(abs(temp(n,:)).^2)^0.5;
                %    temp1(n,:)      = temp1(n,:) / sum(abs(temp1(n,:)).^2)^0.5;
                % end
                temp = Mxy(:,:,end).';
                temp1 = Mxy0(:,:,end).';
                img_sim_fw = temp;
                Sfw_mgre = temp1;
    
                TE_all = [param.TE_flash param.TE_flash param.TE_flash param.TE_flash param.TEs(:)'];
                Npt = param.ncontrast;
                TE_all = reshape(TE_all, 1, Npt);  
                ppm2Hz = single(this.B0*this.gyro);
                phi_dbkg = (2*pi) .* ((ppm2Hz .* freqBKG(:)) * TE_all);
    
                sHat = (img_sim_mw+img_sim_aew+img_sim_fw).*exp(1i.*pini(:)).*exp(1i.*phi_dbkg);%.*exp(1i.*2.*pi.*((freqBKG).').*this.B0.*this.gyro);
                Sreal = real(sHat);
                Simag = imag(sHat);

            else
                % load model
                switch fitting.model
                   case 'invivo_1mm'
                        model_dir = fullfile('/autofs/cluster/berkin/yuting/MATLAB/demo/gacelle-main/ANN_EPGXgen20240927/ablation2/MIMOSA_featAbl_loss1_tau0p01_v1/A05_freqRaw__PT-prepOH__FW0_F16/');   % <-- CHANGE
                        model_fn  = fullfile(model_dir, 'A05_freqRaw__PT-prepOH__FW0_F16_final.mat');
                    case 'invivo_700um'
                        model_dir = fullfile('/autofs/cluster/berkin/yuting/MATLAB/demo/gacelle-main/ANN_EPGXgen20240927/ablation700um/MIMOSA_featAbl_ARCHxFEAT_loss1_tau0p01_v1/A04_NDIV2__freqRaw__PT-prepOH_F16_seed0/');   % <-- CHANGE
                        model_fn  = fullfile(model_dir, 'A04_NDIV2__freqRaw__PT-prepOH_F16_seed0_final.mat'); 
                    case 'sim'
                        model_dir = fullfile('/autofs/cluster/berkin/yuting/MATLAB/demo/gacelle-main/ANN_EPGXgen20240927/ablation2/MIMOSA_sim_featAbl_loss1_tau0p01_v1/A05_freqRaw__PT+prepOH__FW0_F16/');   % <-- CHANGE
                        model_fn  = fullfile(model_dir, 'A05_freqRaw__PT+prepOH__FW0_F16_model.mat');
                    case 'L1'
                        model_dir = fullfile('/autofs/cluster/berkin/yuting/MATLAB/demo/gacelle-main/ANN_EPGXgen20240927/ablation6_loss_stable/MIMOSA_loss8Abl_minChange_v1/L07_LOSS07__RL1_F16_seed1/');   % <-- CHANGE
                        model_fn  = fullfile(model_dir, 'L07_LOSS07__RL1_F16_seed1_final.mat');
                    case 'WRL1'         
                        model_dir = fullfile('/autofs/cluster/berkin/yuting/MATLAB/demo/gacelle-main/ANN_EPGXgen20240927/ablation6_loss_stable/MIMOSA_loss8Abl_minChange_v1/L05_LOSS05__WRL1_F16_seed1/');   % <-- CHANGE
                        model_fn  = fullfile(model_dir, 'L05_LOSS05__WRL1_F16_seed1_final.mat');
                    case 'invivo_IR750V' % for CimaX
                        model_dir = fullfile('/autofs/cluster/berkin/yuting/MATLAB/demo/gacelle-main/ANN_EPGXgen20240927/ablation2/MIMOSA_IR750/A05_freqRaw__PT+prepOH__FW0_F16');
                        model_fn  = fullfile(model_dir, 'A05_freqRaw__PT+prepOH__FW0_F16_model.mat');
                 
                end
                S = load(model_fn, 'dlnet');
                load(model_fn, 'feature_idx');

                dlnet = S.dlnet;
                alpha_leaky = dlnet.alpha;
                parameters  = dlnet.parameters;
                Npt = param.ncontrast;

                param = build_MIMOSA_acq_times(param, param.echo2use);
                param.tmax = max(param.t_acq14_s);

                executionEnvironment = "gpu";
                nfeatures  = length(feature_idx);
                numOutputs = 2;

                % ---- Ranges used in training (must match exactly) ----
                idx = find(ismember(this.modelParams,'T1MW'));
                T1_MW_range_s   = [this.lb(idx), this.ub(idx)];
                idx = find(ismember(this.modelParams,'T1IEW'));
                T1_IEW_range_s  =  [this.lb(idx), this.ub(idx)];
                idx = find(ismember(this.modelParams,'r2MW'));
                R2_MW_range_Hz  = [this.lb(idx), this.ub(idx)];
                idx = find(ismember(this.modelParams,'r2IEW'));
                R2_IEW_range_Hz = [this.lb(idx), this.ub(idx)];

                % TE list (1x14)
                TE_all = [param.TE_flash param.TE_flash param.TE_flash param.TE_flash param.TEs(:)'];
                
                TE_all = reshape(TE_all, 1, Npt);                     % [1, Npt]
                % TE_all_dl = dlarray(TE_all, 'CB');                    % treat as constant; label not super critical
                
                % --- Ensure scalars/maps are dlarray or derived from dlarray ---
                mwf_2d = repmat(mwf(:), 1, Npt);                      % dlarray if mwf is dlarray
                iew_2d = repmat(iew(:), 1, Npt);
                fw_2d  = repmat(fw(:),  1, Npt);
                b1_2d  = repmat(b1_val(:), 1, Npt);
                
                ptIndexNorm = single((0:Npt-1)/(Npt-1));
                pt_2d = repmat(ptIndexNorm, num_vox, 1);              % numeric const ok; can wrap dlarray if you want
                
                % --- Decay features (stay dlarray) ---
                T2s_MW_decay  = exp(-mwf_2d*0 + TE_all .* (-r2sMW(:)));    % keep dlarray path
                T2s_IEW_decay = exp(-iew_2d*0 + TE_all .* (-r2sIEW(:)));
                
                ppm2Hz = single(this.B0*this.gyro);
                phi_mw  = 2*pi * (ppm2Hz .* (freqMW(:))  * TE_all);
                phi_iew = 2*pi * (ppm2Hz .* (freqIEW(:)) * TE_all);
                phi_fw = zeros(num_vox, Npt, 'like', prototype);% FW freq
                phi_dbkg  = 2*pi * (ppm2Hz .* (freqBKG(:))  * TE_all);

                sinPhi_MW  = sin(phi_mw);   cosPhi_MW  = cos(phi_mw);
                sinPhi_IEW = sin(phi_iew);  cosPhi_IEW = cos(phi_iew);
                sinPhi_FW = sin(phi_fw);  cosPhi_FW = cos(phi_fw);

                 % ---- Prep one-hot (4) ----
                prepOH = zeros(num_vox, Npt, 4, 'like', mwf_2d);
                prepOH(:,1,1)=1;
                prepOH(:,2,2)=1;
                prepOH(:,3,3)=1;
                prepOH(:,4,4)=1;

                  % ---- T1/R2 normalized (match training) ----
                T1MWn = (pars.T1MW(:)  - min(T1_MW_range_s))  ./ diff(T1_MW_range_s);
                T1IEWn= (pars.T1IEW(:) - min(T1_IEW_range_s)) ./ diff(T1_IEW_range_s);
                R2MWn = (pars.r2MW(:)   - min(R2_MW_range_Hz)) ./ diff(R2_MW_range_Hz);   % <-- requires pars.r2MW to exist
                R2IEWn= (pars.r2IEW(:)  - min(R2_IEW_range_Hz))./ diff(R2_IEW_range_Hz);
            
                T1MWn_2d  = repmat(T1MWn,  1, Npt);
                T1IEWn_2d = repmat(T1IEWn, 1, Npt);
                R2MWn_2d  = repmat(R2MWn,  1, Npt);
                R2IEWn_2d = repmat(R2IEWn, 1, Npt);

                % ---- T2prep features (40ms,80ms) from T2=1/R2 ----
                T2_MW_s  = 1 ./ (pars.r2MW(:)  + eps('single'));   % <-- requires pars.r2MW
                T2_IEW_s = 1 ./ (pars.r2IEW(:) + eps('single'));
            
                TE_T2prep_s = single([40e-3, 80e-3]);
                T2pMW40 = exp(-TE_T2prep_s(1) ./ (T2_MW_s  + eps('single')));
                T2pMW80 = exp(-TE_T2prep_s(2) ./ (T2_MW_s  + eps('single')));
                T2pIEW40= exp(-TE_T2prep_s(1) ./ (T2_IEW_s + eps('single')));
                T2pIEW80= exp(-TE_T2prep_s(2) ./ (T2_IEW_s + eps('single')));
            
                T2pMW40_2d  = repmat(T2pMW40,  1, Npt);
                T2pMW80_2d  = repmat(T2pMW80,  1, Npt);
                T2pIEW40_2d = repmat(T2pIEW40, 1, Npt);
                T2pIEW80_2d = repmat(T2pIEW80, 1, Npt);

                % ---- T1basisAcq using cached acquisition times ----
                t_acq14_s = single(param.t_acq14_s(:).');      % [1,14]
                t_acq_2d  = repmat(t_acq14_s, num_vox, 1);
            
                T1basisAcq_MW  = exp(-t_acq_2d ./ (pars.T1MW(:)  + eps('single')));
                T1basisAcq_IEW = exp(-t_acq_2d ./ (pars.T1IEW(:) + eps('single')));
            
                t_acq_norm = t_acq_2d ./ (single(param.tmax) + eps('single'));


                % --- Build features by concatenating dlarrays (NO single() casting here) ---
                features = cat(3, ...
                    mwf_2d, ...
                    iew_2d, ...
                    T2s_MW_decay, ...
                    T2s_IEW_decay, ...
                    repmat(freqMW(:),  1, Npt), ...
                    repmat(freqIEW(:), 1, Npt), ...
                    b1_2d, ...
                    pt_2d, ...
                    sinPhi_MW,  cosPhi_MW, ...
                    sinPhi_IEW, cosPhi_IEW, ...
                    prepOH,...
                    T1MWn_2d, ...               
                    T1IEWn_2d, ...              
                    R2MWn_2d, ...              
                    R2IEWn_2d, ...              %
                    T2pMW40_2d, ...            
                    T2pMW80_2d, ...             
                    T2pIEW40_2d, ...            
                    T2pIEW80_2d, ...            %
                    T1basisAcq_MW, ...          %
                    T1basisAcq_IEW, ...         
                    t_acq_norm ...              
                    );
                featuresFull(:,:,1) = mwf_2d;
                featuresFull(:,:,2) = iew_2d;
                featuresFull(:,:,3) = single(T2s_MW_decay);
                featuresFull(:,:,4) = single(T2s_IEW_decay);
                featuresFull(:,:,5) = T1MWn_2d;
                featuresFull(:,:,6) = T1IEWn_2d;
                featuresFull(:,:,7) = R2MWn_2d;
                featuresFull(:,:,8) = R2IEWn_2d;
                featuresFull(:,:,9) = b1_2d;
                
                % 10..11 raw frequency
                featuresFull(:,:,10) = repmat(freqMW(:),  1, Npt);
                featuresFull(:,:,11) =  repmat(freqIEW(:), 1, Npt);
                
                % 12..15 phase trig
                featuresFull(:,:,12) = single(sinPhi_MW);
                featuresFull(:,:,13) = single(cosPhi_MW);
                featuresFull(:,:,14) = single(sinPhi_IEW);
                featuresFull(:,:,15) = single(cosPhi_IEW);
                
                % 16..21 time/context groups
                featuresFull(:,:,16) = pt_2d;

                featuresFull(:,:,17:20) = prepOH;
                
                featuresFull(:,:,21) = t_acq_norm;

                % 22 optional FW
                featuresFull(:,:,22) = fw_2d;
                
                 % ---------------- Subset and output ----------------
                features = featuresFull(:,:,feature_idx);
                
                
                % ---- Forward network ----
                dlin3 = permute(features, [3 1 2]);                 % [31,B,14]
                dlin  = reshape(dlin3, [nfeatures, num_vox*Npt]);   % [31,B*14]
                dlX   = dlarray(dlin, 'CB');
                dlX   = gpuArray(dlX);                                % if you want


                U = mlp_model_leakyRelu(parameters, dlX, alpha_leaky);  % [2, N*14]
                % U = mlp_model_leakyRelu_chunked(parameters, dlX, alpha_leaky);  % [2, N*14]

                S_pred = reshape(U, [2, num_vox, Npt]);                 % [2, B, 14]

                Re0 = squeeze(S_pred(1,:,:));                            % [B,14]
                Im0 = squeeze(S_pred(2,:,:));                            % [B,14]


                cb = cos(phi_dbkg);  sb = sin(phi_dbkg); 
                % implicit expansion to [B,14]
                Re = Re0 .* cb - Im0 .* sb;
                Im = Re0 .* sb + Im0 .* cb;

                
                % scale by S0 (still real)
                Re = Re .* S0(:);
                Im = Im .* S0(:);
                
                % apply phase offset pini using real rotation (NO complex numbers)
                p = pini(:);                                            % [B,1] (should be real)
                cp = cos(p);  sp = sin(p);                              % [B,1]

                
                % implicit expansion to [B,14]
                Sreal = Re .* cp - Im .* sp;
                Simag = Re .* sp + Im .* cp;
                
            end

            if ~fitting.isComplex
                s = sqrt(Sreal.^2 + Simag.^2);
            else
                s = cat(5,Sreal,Simag);
            end

            % vectorise to match maksed measurement data
            % s = utils.reshape_ND2GD(s,[]);
            s = utils.reshape_ND2GD(permute(s,[1,3,4,2,5]),[]);% yc: the format of s is voxel,t,1,1,2

        end
        
        %% Utilities

        % validate extra data
        function [extraData,mask] = validate_data(this,data,extraData,mask,fitting)

           % % check if the signal is monotonic decay
           %  [~,I] = max(abs(data),[],4);
           %  mask = and(mask,I<4);

            dims = size(data,1:3);

            if ~fitting.DIMWI.isFitIEW && ~isfield(extraData,'IEW')
                error('Field IEW is missing in exraData structure variable for DIMWI model');
            end
            if ~isfield(extraData,'freqBKG')
                extraData.freqBKG = zeros(dims);
                if fitting.isComplex
                    warning('No total field map is provided for fitting complex-valued data.');
                end
            end
            if ~isfield(extraData,'pini')
                % extraData.pini = zeros(dims);
                extraData.pini = angle( data(:,:,:,1) ./ exp(1i* 2*pi*extraData.freqBKG * (this.B0*this.gyro) .* permute(this.te(1),[2 3 4 1])));
            end

            fields = fieldnames(extraData); for kfield = 1:numel(fields); extraData.(fields{kfield}) = single( extraData.(fields{kfield})); end
            
            % thresholding based on single compartment R2*
            % [R2s0,~]    = this.R2star_trapezoidal(abs(data(5:end,:)),this.te);% yc, use for GD only last mgre part
            [R2s0,~]    = this.R2star_trapezoidal(abs(data(:,:,:,5:end)),this.te);% yc, use for ND only last mgre part
           
            mask        = and(mask,R2s0>this.thres_R2star);

            % DIMWI
            if ~fitting.DIMWI.isFitFreqMW || ~fitting.DIMWI.isFitFreqIEW || ~fitting.DIMWI.isFitR2sIEW
                % fibre fraction
                if isfield(extraData,'ff')
                    extraData.ff                        = bsxfun(@rdivide,extraData.ff,sum(extraData.ff,4));
                    mask                                = and(mask,min(~isnan(extraData.ff),[],4));
                    extraData.ff(isnan(extraData.ff))   = 0;
                else
                    error('Fibre fraction map is required for DIMWI!');
                end
                % fibre orientation
                if ~isfield(extraData,'theta')
                    if ~isfield(extraData,'fo')
                        error('Fibre orientation map is required for DIMWI!');
                    else
                        fo    = double(extraData.fo); % fibre orientation w.r.t. B0
                        theta = zeros(size(extraData.ff));
                        for kfo = 1:size(fo,5)
                            theta(:,:,:,kfo) = this.AngleBetweenV1MapAndB0(fo(:,:,:,:,kfo),this.B0dir);
                        end
                        extraData.theta = single(theta);
                        extraData = rmfield(extraData,"fo");
                    end
                end
            else
                extraData.theta = zeros(dims,'single');
                extraData.ff    = ones(dims,'single');
            end

        end

        % normalise input data based on masked signal intensity at 98%
        function [img, scaleFactor] = prepare_data(this,img,extraData, mask) % yc: revise to use PD as norm factor

            % [~,S0] = this.R2star_trapezoidal(abs(img(5:end,:)),this.te);%
             % [~,S0] = this.R2star_trapezoidal(abs(img(:,:,:,5:end)),this.te);% yc, for ND format
            S0 = extraData.PD;
            scaleFactor = prctile( S0(mask), 98);

            img = img ./ scaleFactor;

        end
        
    end

    methods(Static)

        % compute weights for optimisation
        function w = compute_optimisation_weights(data,fitting)
        % 
        % Output
        % ------
        % w         : ND signal masked wegiths that matches the arrangement in masked data later on
        %
            if fitting.isWeighted
                switch lower(fitting.weightMethod)
                    case 'norm'
                       % weights using echo intensity, as suggested in Nam's paper
                        w = sqrt(abs(data));
                    case '1stecho'
                        p = fitting.weightPower;
                        % weights using the 1st echo intensity of each flip angle
                        w = bsxfun(@rdivide,abs(data).^p,abs(data(:,:,:,1)).^p);
                end
            else
                w = ones(size(data));
            end

            w(w>1) = 1; w(w<0) = 0;
            
            % separate real/imaginary parts into 6th dim
            if fitting.isComplex
                w = repmat(w,1,1,1,1,2);
            end
        end
       
        %% signal
        % simple 2-pool matrix inversion
        function [m0,mwf] = superfast_mwi_2m_standard(img,te,t2s)
        %
        % Input
        % --------------
        % img           : multi-echo GRE image, 4D [row,col,slice,TE]
        % te            : echo times in second
        % t2s           : T2* of the two pools, in second, [T2sMW,T2sIEW], if empty
        %                 then default values for 3T will be used
        %
        % Output
        % --------------
        % m0            : proton density of each pool, 4D [row,col,slice,pool]
        % mwf           : myelin water fraction map, range [0,1]
        %
        % Description:  Direct matrix inversion based on simple 2-pool model, i.e.
        %               S(te) = E2s * M0
        %               Useful to estimate initial starting points for MWI fitting
        %
        % Kwok-shing Chan @ DCCN
        % k.chan@donders.ru.nl
        % Date created: 13 Nov 2020
        % Date modified:
        %
        %

            % get size in all image dimensions
            dims = size(img,1:3);
            
            % check input
            if isempty(t2s)
                t2s = [10e-3, 60e-3];   % 3T, [MW, IEW], in second
            end
            
            % T2* decay matrix
            E2s1    = exp(-te(:)/t2s(1));
            E2s2	= exp(-te(:)/t2s(2));
            E2s     = [E2s1,E2s2];
            
            tmp = reshape(abs(img),prod(dims),length(te));
            
            m0 = E2s \ tmp.';
            m0 = reshape(m0.',[dims length(t2s)]);
            
            % compute MWF
            mwf = m0(:,:,:,1) ./ sum(m0,4);
            mwf(mwf<0)      = 0;
            mwf(mwf>1)      = 1;
            mwf(isnan(mwf)) = 0;
            mwf(isinf(mwf)) = 0;
            
            m0(m0 < 0)      = 0;
            m0(isinf(m0))   = 0;
            m0(isnan(m0))   = 0;
        
        end

        % closed form single compartment solution
        function [R2star,S0] = R2star_trapezoidal(img,te)
            % disgard phase information 

            img = double(abs(img));
            te  = double(te);
            
            dims = size(img);

            % yc, add when input is GD format
            % if size(dims)<3
            %     img = reshape(img.',[],1,1,10);% t x voxel to x,y,z,t
            % end
            % dims = size(img);
            
            % main
            % Trapezoidal approximation of integration
            temp = 0;
            for k=1:dims(4)-1
                temp = temp+0.5*(img(:,:,:,k)+img(:,:,:,k+1))*(te(k+1)-te(k));
            end
            
            % very fast estimation
            t2s = temp./(img(:,:,:,1)-img(:,:,:,end));
                
            R2star = 1./t2s;

            S0 = img(1:(numel(img)/dims(end)))'.*exp(R2star(:)*te(1));
            if numel(S0) ~=1
                S0 = reshape(S0,dims(1:end-1));
            end
        end

        %% Utilities
        % check and set default fitting algorithm parameters
        function fitting2 = check_set_default(fitting,data)
            % get basic fitting setting check
            fitting2 = askadam.check_set_default_basic(fitting);

            % check weighted sum of cost function
            if ~isfield(fitting,'isWeighted');      fitting2.isWeighted     = true;         end
            if ~isfield(fitting,'weightMethod');    fitting2.weightMethod   = '1stecho';    end
            if ~isfield(fitting,'weightPower');     fitting2.weightPower    = 1;            end
            
            % check hollow cylinder fibre model parameters
            if ~isfield(fitting,'DIMWI') || ~isfield(fitting.DIMWI,'isFitFreqMW');  fitting2.DIMWI.isFitFreqMW  = true; end
            if ~isfield(fitting,'DIMWI') || ~isfield(fitting.DIMWI,'isFitFreqIEW');  fitting2.DIMWI.isFitFreqIEW  = true; end
            if ~isfield(fitting,'DIMWI') || ~isfield(fitting.DIMWI,'isFitR2sIEW');   fitting2.DIMWI.isFitR2sIEW   = true; end
            if ~isfield(fitting,'DIMWI') || ~isfield(fitting.DIMWI,'isFitIEW');     fitting2.DIMWI.isFitIEW     = true; end

            % get customised fitting setting check
            if ~isfield(fitting,'regmap');      fitting2.regmap = 'MWF'; end
            if ~isfield(fitting,'start');       fitting2.start  = 'prior'; end

            if ~isfield(fitting,'isComplex');   fitting2.isComplex = true; end
            if isreal(data);                    fitting.isComplex = false;  end

        end

        function display_algorithm_info(fitting)
            %%%%%%%%%% 3. display some algorithm parameters %%%%%%%%%%
            disp('--------------');
            disp('Fitting option');
            disp('--------------');
            % type of fitting
            if fitting.isComplex
                disp('Fitting with complex-valued data');
            else 
                disp('Fitting with magnitude data');
            end

            disp('Cost function options:');
            if fitting.isWeighted
                disp('Cost function weighted by echo intensity: True');
                disp(['Weighting method: ' fitting.weightMethod]);
                if strcmpi(fitting.weightMethod,'1stEcho')
                    disp(['Weighting power: ' num2str(fitting.weightPower)]);
                end
            else
                disp('Cost function weighted by echo intensity: False');
            end

            disp('------------------------------------');
            disp('Diffusion informed MWI model options');
            disp('------------------------------------');
            if ~fitting.DIMWI.isFitIEW
                disp('Fit intra-axonal volume fraction  : False');
            else
                disp('Fit intra-axonal volume fraction  : True');
            end
            if ~fitting.DIMWI.isFitFreqMW
                disp('Fit frequency - myelin water      : False');
            else
                disp('Fit frequency - myelin water      : True');
            end
            if ~fitting.DIMWI.isFitFreqIEW
                disp('Fit frequency - intra-axonal water: False');
            else
                disp('Fit frequency - intra-axonal water: True');
            end
            if ~fitting.DIMWI.isFitR2sIEW
                disp('Fit R2* - extra-cellular water    : False');
            else
                disp('Fit R2* - extra-cellular water    : True');
            end

            disp('------------------------------------');

        end

        function theta = AngleBetweenV1MapAndB0(v1,b0dir)
        %
        % Input
        % --------------
        % v1            : 4D fibre orientation map in vector form
        % b0dir         : 1D vector of B0 direction
        %
        % Output
        % --------------
        % theta         : 3D angle map, in rad
        %
        % Description:
        %
        % Kwok-shing Chan @ DCCN
        % k.chan@donders.ru.nl
        % Date created: 20 March 2019
        % Date last modified: 25 October 2019
        %
        %

            % replicate B0 direction to all voxels
            b0dirmap = permute(repmat(b0dir(:),1,size(v1,1),size(v1,2),size(v1,3)),[2 3 4 1]);
            % compute angle between B0 direction and fibre orientation
            theta = atan2(vecnorm(cross(v1,b0dirmap),2,4), dot(v1,b0dirmap,4));
            
            % make sure the angle is in range [0, pi/2]
            theta(theta> (pi/2)) = pi - theta(theta> (pi/2));
        
        end
    
    end
end