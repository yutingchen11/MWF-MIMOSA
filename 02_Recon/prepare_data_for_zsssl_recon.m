
%--------------------------------------------------------------------------
%% Read twix data and sort kspace based on seq file
%--------------------------------------------------------------------------

addpath(genpath('utils'));

clear; clc;
flag.save_kspace = 1;
flag.save_nii = 1;

data_file_path='rawdata/meas_MID00348_FID00459_mwf_mimosa_R8.dat';

[p,n,e] = fileparts(data_file_path);
basic_file_path=fullfile(p,n);

twix_obj = mapVBVD(data_file_path);
data_unsorted = twix_obj{end}.image.unsorted();
[adc_len,ncoil,readouts]=size(data_unsorted);


% Read params from seq file
pulseq_file_path = ['../01_Pulse_Seq' '/MWF_MIMOSA_1iso' '.seq'];
seq=mr.Sequence();
seq.read(pulseq_file_path);


N = seq.getDefinition('Matrix');
nTR = seq.getDefinition('nTR');
nETL = seq.getDefinition('nETL');
os_factor=1;
traj_y = seq.getDefinition('traj_y');
traj_z = seq.getDefinition('traj_z');
step_size = 4 * nETL;
%#######
TEs = seq.getDefinition('TES_mte');
nechoes = seq.getDefinition('num_echoes');
esp_mte = seq.getDefinition('esp_mte');
TR_mte = seq.getDefinition('TR_mte');
%###########
nacq = 4+10;
%% prescan
ny_lr = 32;
nz_lr = 32;
k_ref = reshape(permute(data_unsorted(:,:,1:ny_lr*nz_lr),[1,3,2]),N(1),ny_lr,nz_lr,ncoil);
jimg2(rsos(ifft3call(k_ref),4));
save('ref_mimosa_R8_bp.mat','k_ref');
data_unsorted_data = data_unsorted(:,:,ny_lr*nz_lr+1:end);
clear data_unsorted
%% Pre-allocate
img4d_data = zeros([N(1)*os_factor N(2) N(3) 4+nechoes]);
kspace = zeros([N(1)*os_factor N(2) N(3)  ncoil 4+nechoes]);
mask_traj = zeros([N(2) N(3) nETL 3+nechoes]);
%%######
ind14 = [];
indmte = [];
for t = 0:nTR-1
    c1_start =  (4 + nechoes)*t*nETL+1;
    c4_end = (4 + nechoes)*t*nETL + 4*nETL;
    c5_start = c4_end+1;
    c10_end = (4 + nechoes)*(t+1)*nETL;
    ind14 = [ind14, c1_start:c4_end];
    indmte = [indmte, c5_start:c10_end];
end
data_c1_to_c4 = data_unsorted_data(:,:,ind14(:));% extract data of each contrast
data_c5_to_c10 = data_unsorted_data(:,:,indmte(:));% extract data of each contrast
readout_c4 = size(data_c1_to_c4,3);

for contrast = 1:4
    indices = [];
    for segment_start = (1+nETL*(contrast-1)):step_size:readout_c4
        segment_end = min(segment_start + nETL - 1, readout_c4);
        indices = [indices, segment_start:segment_end];
    end

    data = data_c1_to_c4(:,:,indices(:));% extract data of each contrast
    kspace_contrast = zeros([N(1)*os_factor N(2) N(3) ncoil]);

    ii = 1;
    for rr= 1:nTR
        for ee = 1:nETL
            ky = traj_y((contrast-1)*nTR*nETL+ee+(rr-1)*nETL);
            kz = traj_z((contrast-1)*nTR*nETL+ee+(rr-1)*nETL);
%             disp((contrast-1)*nTR*nETL+ee+(rr-1)*nETL);
            if sum(kspace_contrast(:,ky,kz,:),'all')==0
                kspace_contrast(:,ky,kz,:) = data(:,:,ii);% To avoid overwriting
                mask_traj(ky,kz,ee,contrast) = 1;
            else
                disp(['Skip ky=' num2str(ky) ',kz=' num2str(kz)])
            end
            ii = ii+1;
             % figure(contrast);hold on;scatter(ky,kz,10,contrast);
             
        end
    end

    im = fft3c2(kspace_contrast);
    im3D = abs(sum(im.*conj(im),ndims(im))).^(1/2);
    % jimg2(im3D);
    img4d_data(:,:,:,contrast) = im3D;
    kspace(:,:,:,:,contrast) = kspace_contrast;
end

%%
for contrast = 5:nechoes+4

    ind_acq5 = [];
    for i = 1:nTR
        ind_acq5 = [ind_acq5, (contrast-4)+nETL*nechoes*(i-1):nechoes:nETL*nechoes*i];
    end

    data = data_c5_to_c10(:,:,ind_acq5(:));
    kspace_contrast = zeros([N(1)*os_factor N(2) N(3) ncoil]);
    % for no cplm, outc
         ii = 1;
    for nt = 1:nTR
        for yy=1:nETL
            ind_traj = 5;
            index = (ind_traj-1)*nTR*nETL+(nt-1)*nETL+yy;
            index_oc =  (ind_traj-1)*nTR*nETL+(nt)*nETL - yy +1;
            ky = traj_y(index_oc);
            kz = traj_z(index_oc);
            if sum(kspace_contrast(:,ky,kz,:),'all')==0
                if mod(contrast,2)==0
                    kspace_contrast(:,ky,kz,:) = data(end:-1:1,:,ii);% To avoid overwriting
                else
                    kspace_contrast(:,ky,kz,:) = data(:,:,ii);% To avoid overwriting
                end
                mask_traj(ky,kz,ee,ind_traj) = 1;
                
            else
                disp(['Skip ky=' num2str(ky) ',kz=' num2str(kz)])
            end
            ii = ii+1;
%             mask_traj(ky,kz,contrast) = 1;
%             mask_traj(ky,kz,contrast) = mask_traj(ky,kz,contrast)+1;
        end
    end

    im = fft3c2(kspace_contrast);
    im3D = abs(sum(im.*conj(im),ndims(im))).^(1/2);
    % jimg2(im3D)

    img4d_data(:,:,:,contrast) = im3D;
    kspace(:,:,:,:,contrast) = kspace_contrast;
end
% jimg2(im3D)

%%
clearvars -except kspace data_file_path

size_data = size(kspace(:,:,:,1,1));

%--------------------------------------------------------------------------
%% patref scan
%--------------------------------------------------------------------------
   
load('ref_mimosa_R8_bp.mat');
ref = flip(flip(k_ref, 2), 3);

img_ref = ifft3call(ref);

imagesc3d2(rsos(img_ref,4), s(img_ref)/2, 1, [0,0,0], [-0,3e-3]), setGcf(.5)


%--------------------------------------------------------------------------
%% coil compression
%--------------------------------------------------------------------------

num_chan = 20;  % num channels to compress to

[ref_svd, cmp_mtx] = svd_compress3d(ref, num_chan, 1);  

rmse(rsos(ref_svd,4), rsos(ref,4))

N = size(kspace(:,:,:,1,1));
num_eco = size(kspace,5);

kspace_svd = zeross([N,num_chan,num_eco]);

for t = 1:size(kspace,5)
    kspace_svd(:,:,:,:,t) = svd_apply3d(kspace(:,:,:,:,t), cmp_mtx);
end

rmse(rsos(kspace_svd,4), rsos(kspace,4))
clear kspace

%--------------------------------------------------------------------------
%% interpolate patref by zero padding to the high res matrix size
%--------------------------------------------------------------------------

size_data = size(kspace_svd(:,:,:,1,1));
size_patref = size(ref_svd(:,:,:,1,1));
% pad ref to be the same size
patref_pad = padarray( ref_svd, [size_data-size_patref, 0, 0, 0]/2 );

img_patref_pad = ifft3c(patref_pad);

imagesc3d2( rsos(img_patref_pad,4), s(img_patref_pad)/2, 10, [0,0,0], [0,2e-4])
%--------------------------------------------------------------------------
%% calculate sens map using ESPIRiT: parfor
%--------------------------------------------------------------------------
num_acs = min(size_patref);
kernel_size = [6,6];
eigen_thresh = 0.7;

receive = zeross(size(kspace_svd(:,:,:,:,1)));


delete(gcp('nocreate'))
c = parcluster('local');    

total_cores = c.NumWorkers;  
parpool(ceil(total_cores/4))

    
tic
parfor slc_select = 1:s(img_patref_pad,1)     
    disp(num2str(slc_select))
    
    [maps, weights] = ecalib_soft( fft2c( sq(img_patref_pad(slc_select,:,:,:)) ), num_acs, kernel_size, eigen_thresh );

    receive(slc_select,:,:,:) = permute(dot_mult(maps, weights >= eigen_thresh ), [1,2,4,3]);
end 
toc
receive = abs(receive) .* exp(1i * angle( receive .* repmat(conj(receive(:,:,:,1)), [1,1,1,num_chan]) ));

 
save([pwd, '/receive_svd_', num2str(num_chan), 'ch_mimosa_R8_bp.mat'], 'receive', '-v7.3')


%% SENSE
%multi slice
 % addpath(genpath('F:\beifen\MGH_new\MATLAB\demo\JLORAKS_2D-4Yuting'))
mask_all = kspace_svd ~= 0;
[nx, ny, nz, nc, ne] = size(kspace_svd);
% sense = zeros([nx,ny,nz,ne]);

slice_idx = 150; %1:nx
    disp(['slice:' num2str(slice_idx)]);
    tmp = fftshift(ifft(ifftshift(kspace_svd),[],1))* sqrt(nx);   % send kx into image domain
    kData = squeeze(tmp(slice_idx,:,:,:,:));        % pick single slice in x direction
    
    kMask = squeeze(mask_all(slice_idx,:,:,:,:));
    
    coil_sens = single(squeeze(receive(slice_idx,:,:,:)));
    %%%%%%%%%%%%%%%%
    
    data = kMask .* kData;



    fig_num = 3;    % figure number

    Ah = @(x) vect(sum( conj(coil_sens) .* ift2(  kMask .* reshape(x, [ny nz nc ne]) ), 3  ));
    AhA = @(x)   vect( sum( conj(coil_sens).*ift2( kMask.*ft2(coil_sens.* reshape(x, [ny nz 1 ne]))),3));


    Ahd = Ah(data);


    disp('SENSE Reconstruction');
    tic
    [recon,flag,relres,iter] = pcg( AhA, Ahd, [], 10); % iters
    toc
    recon = reshape(recon,[ny nz ne]);
    sense = recon;


mosaic(sq(sense), 2, 5, fig_num, '', [0, 2e-4],180)

%% generate data for reconstruction
mask_all = kspace_svd ~= 0;
[nx, ny, nz, nc, ne] = size(kspace_svd);

tmp = fftshift(ifft(ifftshift(kspace_svd),[],1))* sqrt(nx); 
%% slice by slice
if ~exist('zsssl_recon/data', 'dir')
    mkdir('zsssl_recon/data'); 
end
cd('zsssl_recon/data/')

for ss = 50:230 % neck was skipped here, use 1:end if needed
    kspace =  squeeze(tmp(ss,:,:,:,:));
    sens_maps = squeeze(receive(ss,:,:,:));
    mask = sq(mask_all(:,:,1,:));
    filename = sprintf('mimosa_R8_slc_%03d', ss);
    save([filename '.mat'],'kspace','sens_maps','mask')
end
