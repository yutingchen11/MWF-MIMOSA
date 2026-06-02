% 1. generate dict first by run gen_MWF_MIMOSA_dict_1mm.m
% 2. install sepia Tool: https://github.com/kschan0214/sepia
% 3. install χ-separation Tool: https://github.com/SNU-LIST/chi-separation
% including 
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

% Set x-separation tool directory path
%%
addpath(genpath('/autofs/cluster/berkin/yuting/MATLAB/demo/chi-separation-main/xiangruili-dicm2nii-3fe1a27'))
% Tools for NIfTI and ANALYZE image (https://kr.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
addpath(genpath('/autofs/cluster/berkin/yuting/MATLAB/demo/chi-separation-main/NIfTI_20140122'))

% Download onnxconverter Add-on, and then install it.
% Deep Learning Toolbox Converter for ONNX Model Format 
% (https://kr.mathworks.com/matlabcentral/fileexchange/67296-deep-learning-toolbox-converter-for-onnx-model-format)

% Set QSM tool directory path 
% STI Suite (Version 3.0) (https://people.eecs.berkeley.edu/~chunlei.liu/software.html)
addpath(genpath('/autofs/cluster/berkin/yuting/MATLAB/demo/chi-separation-main/STISuite_V3.0/STISuite_V3.0'))

% MEDI toolbox (http://pre.weill.cornell.edu/mri/pages/qsm.html)
addpath(genpath('/autofs/cluster/berkin/yuting/MATLAB/demo/chi-separation-main/MEDI_toolbox'))

% SEGUE toolbox (https://xip.uclb.com/product/SEGUE)
addpath(genpath('/autofs/cluster/berkin/yuting/MATLAB/demo/chi-separation-main/SEGUE_28012021'))

% mritools toolbox (https://github.com/korbinian90/CompileMRI.jl/releases)
addpath(genpath('/autofs/cluster/berkin/yuting/MATLAB/demo/mritools_ubuntu-20.04_4.5.3'))

% Set x-separation tool directory path
home_directory = '/autofs/cluster/berkin/jiye/ProcessingCode/Chisep_Toolbox_v1.2.1_pfile';
addpath(genpath(home_directory))
%%
addpath(genpath('gacelle-mwf-mimosa'));  
addpath(genpath('utils'));  
%% load dict
dict = h5read('dict/dict_mwf_mimosa_bp_50b1_126T2s_double.h5','/dict');
t2s_entries  = [1:1:100, 105:5:200 210:50:500];
%% load look up table
load('dict/ielookup_4mimosa.mat');
%% load recon data
load('../02_Recon/zsssl_recon_1mm/recon/mwf_mimosa_R8_zsssl_tfv2.mat');

path_B1 = '../02_Recon/rawdata/dcm_b1_tfl_b1map_20260302130524_6001.nii';
%--------------------------------------------------------------------------
%% if save nifti files
is_save_nii = 1;
%--------------------------------------------------------------------------
%% set protocol for 1mm iso
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


img = abs(img_all);

%--------------------------------------------------------------------------
%% Load b1 map and resize
%--------------------------------------------------------------------------

matrix_size = [size(img,1),size(img,2),size(img,3)];

img_b1_load = niftiread(path_B1);
img_b1 = double(img_b1_load)/800;       % reference fa 800
% img_b1 = flip(img_b1, 1);  
img_b1 = permute(img_b1, [2 1 3]); 
img_b1 = flip(img_b1,2); %  

img_b1 = imresize3(img_b1, [240,240,195]); %  to match img size, use the one in your protocol
img_b1 = img_b1(:,9:end-8,2:end-2);

imagesc3d2(img_b1, s(img_b1)/2, 51, [-0,-0,0], [0,2])
imagesc3d2(img(:,:,:,1), s(img(:,:,:,1))/2, 52, [-0,-0,0], [0,2e-4])

save('B1_4mimosa.mat','img_b1')
%--------------------------------------------------------------------------
%% single compartment dict fitting for T1/T2/T2*/PD/IE mapping
%--------------------------------------------------------------------------


single_comp_dict_fitting;


imagesc3d2(sq(T1_map), s(sq(T1_map(:,:,:)))/2+[30 0 0], 1, [180 180 180], [0 3000]);colormap hot
imagesc3d2(sq(T2_map), s(sq(T1_map(:,:,:)))/2+[30 0 0], 2, [180 180 180], [0 200]);colormap hot
imagesc3d2(sq(T2s_map), s(sq(T1_map(:,:,:)))/2+[30 0 0], 3, [180 180 180], [0 200]);colormap hot
imagesc3d2(sq(PD_map), s(sq(T1_map(:,:,:)))/2+[30 0 10], 4, [180 180 180], [0 max(abs(PD_map(:)))]);colormap gray
imagesc3d2(sq(IE_map), s(sq(T1_map(:,:,:)))/2+[30 0 10], 5, [180 180 180], [0 1]);colormap gray


save('mapping_mwf_mimosa_R8_bp_zsssl.mat','T1_map','T2_map','T2s_map','PD_map','IE_map','img_all','-v7.3')
% save nii files
if is_save_nii
    nii_T1 = make_nii(sq((T1_map)));
    save_nii(nii_T1,'mimosa_R8_T1.nii');
    nii_T2 = make_nii(sq((T2_map)));
    save_nii(nii_T2,'mimosa_R8_T2.nii');
    nii_T2s = make_nii(sq((T2s_map)));
    save_nii(nii_T2s,'mimosa_R8_T2s.nii');
    nii_pd = make_nii(sq((PD_map)));
    save_nii(nii_pd,'mimosa_R8_pd.nii');
end

%--------------------------------------------------------------------------
%% GACELLE based MWF mapping
%--------------------------------------------------------------------------
MWF_init_optm;% initialization

single_slc_check = 1;% single slc quick check 

MWF_est_GACELLE;

save('gacelle_sub5_A05.mat','results_3D','-v7.3')

if is_save_nii
    MWF = single(sq( results_3D.MWF));
    nii = make_nii(MWF);   
    nii.hdr.dime.datatype = 16;          
    nii.hdr.dime.bitpix = 32;
    save_nii(nii,'mimosa_R8_MWF.nii')
    
    IEW = single(sq(results_3D.IEW));
    nii = make_nii(IEW);   
    nii.hdr.dime.datatype = 16;          
    nii.hdr.dime.bitpix = 32;
    save_nii(nii,'mimosa_R8_IEW.nii')
    
    FW = single(results_3D.FW);
    nii = make_nii(FW);   
    nii.hdr.dime.datatype = 16;         
    nii.hdr.dime.bitpix = 32;
    save_nii(nii,'mimosa_R8_FW.nii')
end

imagesc3d2( MWF, s(MWF)/2+[20 0 0], 12, [180,180,180], [0 0.3])
imagesc3d2( IEW, s(IEW)/2+[20 0 0], 13, [180,180,180], [0 1])
imagesc3d2( FW, s(FW)/2+[20 0 0], 14, [180,180,180], [0 1])


%--------------------------------------------------------------------------
%% chi separation
%--------------------------------------------------------------------------

chi_sep;

imagesc3d2( x_para, s(x_para)/2+[20 0 0], 22, [180,180,180], [0 0.1])
imagesc3d2( x_dia, s(x_dia)/2+[20 0 0], 23, [180,180,180], [0 0.1])
imagesc3d2( x_tot, s(x_tot)/2+[20 0 0], 24, [180,180,180], [-0.1 0.1])

save(['QSM_mwf_mimosa_R8_zsssl_sepnet.mat'],'x_para','x_dia','x_tot')

if is_save_nii
    nii_T1 = make_nii(sq((x_para)));
    save_nii(nii_T1,'mimosa_R8_x_para.nii');
    nii_T2 = make_nii(sq((x_dia)));
    save_nii(nii_T2,'mimosa__R8_x_dia.nii');
    nii_pd = make_nii(sq((x_tot)));
    save_nii(nii_pd,'mimosa_R8_x_tot.nii');
end