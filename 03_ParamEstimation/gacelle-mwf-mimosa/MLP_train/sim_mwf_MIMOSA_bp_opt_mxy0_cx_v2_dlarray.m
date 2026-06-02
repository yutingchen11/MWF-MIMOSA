function [Mz_mtx, Mxy_mtx, Mxy_mgre0] = sim_mwf_MIMOSA_bp_opt_mxy0_cx_v2_dlarray(TR, alpha_deg, esp, turbo_factor, t1_vals, t2_vals, num_reps, echo2use, TR_mte,esp_mte,TEs, t2s_vals, gap_between_readouts, time2relax_at_the_end, b1, inv_eff,local_field_sim_Hz,pd_map )
%
% Modified for dlarray support (Automatic Differentiation)
% v2: considering T2* and phase decay during FLASH readout


prototype_var = t1_vals; 
% -----------------------------

if exist('pd_map', 'var')
    num_voxels = length(pd_map);
else
    num_voxels =  length(t1_vals);
end

if nargin < 18
    M0 = ones(num_voxels,1, 'like', prototype_var);
else
    M0 = pd_map; 
end


if nargin < 16
    inv_eff = 1;
end
if nargin < 15
    b1 = 1;
end
if nargin < 14
    time2relax_at_the_end = 0;
end
if nargin < 13
    gap_between_readouts = 900e-3;
end


relax_T2 = @(Mz, TE_T2prep, T2)  Mz .* exp(-TE_T2prep ./ T2);
inverse_eff = @(M0, IE)  -1.*M0 .* IE;
relax_T2_b1 = @(Mz, TE_T2prep, T2, b1, T1)  Mz .* ( sin(b1 .* pi/2).^2 .* exp(-TE_T2prep ./ T2) + cos(b1 .* pi/2).^2 .* exp(-TE_T2prep ./ T1) );
relax_T1 = @(M0, Mz, delta_t, T1)  M0 - (M0 - Mz) .* exp(-delta_t ./ T1);
relax_T2s = @(Mxy, TE1, T2s)  Mxy .* exp(-TE1 ./ T2s);

relax_T2s_phase_evo = @(Mxy, TE1, T2s, del_f)  Mxy .* exp(-TE1 ./ T2s) .* exp(1i.*2.*pi.*del_f.*TE1);

etl = turbo_factor * esp;
tf = turbo_factor;       
alpha = b1 .* alpha_deg * pi/180;

TE_flash = 2.29e-3;

num_acq = 4 + length(TEs);

% --- Timings Setup (Same as original) ---
delT_M1m_M2m = 29.7e-3;              
delT_M0_M1m = gap_between_readouts - etl - delT_M1m_M2m;
delT_M2m_M3m = etl;                
delT_M1_M2 = 89.7e-3;                
delT_M3m_M1 = gap_between_readouts - etl - delT_M1_M2;
delT_M2_M3 = etl;                   
delT_M2_M6 = gap_between_readouts;         
delT_M4_M5 = 12.8e-3;            
delT_M5_M6 = 100e-3 - 6.45e-3;     
delT_M3_M4 = delT_M2_M6 - delT_M2_M3 - delT_M4_M5 - delT_M5_M6;    
delT_M6_M7 = etl;                   
delT_M7_M8 = 1e-3;
delT_M8_M9 = etl;                   
delT_M9_M10 = 1e-3;
delT_M10_M11 = etl;                
delT_M11_M12 = gap_between_readouts - etl;        
delT_M12_M13 = etl;                 
delT_M12_M13_mte = TR_mte*turbo_factor;

total_event_duration = delT_M0_M1m + delT_M1m_M2m + delT_M2m_M3m + delT_M3m_M1 + delT_M1_M2 + delT_M2_M3 + delT_M3_M4 + delT_M4_M5 + delT_M5_M6 + delT_M6_M7 + delT_M7_M8 + delT_M8_M9 + delT_M9_M10  + delT_M12_M13_mte;

% delT_M13_2end = max(TR - total_event_duration, 0);
delT_M13_2end = 53.5e-3;
if time2relax_at_the_end > 0
    delT_M13_2end = delT_M13_2end + time2relax_at_the_end;
end


Mz_mtx = zeros([14, num_voxels, num_reps], 'like', prototype_var);
Mxy_mtx = zeros([num_acq, num_voxels, num_reps], 'like', prototype_var);
Mxy_mgre0 = zeros([1, num_voxels, num_reps], 'like', prototype_var);

Mstart = M0;

for reps = 1:num_reps
    
    %%%%%% add contrast 0
    M1m = relax_T1(M0, Mstart, delT_M0_M1m, t1_vals);
    M2m = relax_T2_b1(M1m, delT_M1m_M2m - 9.7e-3, t2_vals, b1, t1_vals);
     
    % acq0
    Mz = M2m;
    

    Mxy = zeros(tf, num_voxels, 'like', prototype_var);

    time = 9.7e-3;

    for q = 1:tf
        Mz = relax_T1(M0, Mz, time, t1_vals);

        Mxy(q,:) = relax_T2s_phase_evo(sin(alpha(:,1)) .* Mz,TE_flash,t2s_vals,local_field_sim_Hz);
        Mz = cos(alpha(:,1)) .* Mz;
        time = esp;
    end
 
    M3m = Mz;
    Mxy_acq0 = Mxy;

    %%% start of ori qalas
    M1 = relax_T1(M0, M3m, delT_M3m_M1, t1_vals);
        
    M2 = relax_T2_b1(M1, delT_M1_M2 - 9.7e-3, t2_vals, b1, t1_vals);
     
    % acq1
    Mz = M2;
    Mxy = zeros(tf, num_voxels, 'like', prototype_var); % FIX

    time = 9.7e-3;

    for q = 1:tf
        Mz = relax_T1(M0, Mz, time, t1_vals);
        Mxy(q,:) = relax_T2s_phase_evo(sin(alpha(:,2)) .* Mz,TE_flash,t2s_vals,local_field_sim_Hz);
        Mz = cos(alpha(:,2)) .* Mz;
        time = esp;
    end
 
    M3 = Mz;
    Mxy_acq1 = Mxy;
    
    M4 = relax_T1(M0, M3, delT_M3_M4, t1_vals);

    % inversion efficiency (lookup table or scalar)
    M5 = inverse_eff(M4,inv_eff);

    M6 = relax_T1(M0, M5, delT_M5_M6, t1_vals);
    
    % acq2
    Mz = M6;
    Mxy = zeros(tf, num_voxels, 'like', prototype_var); % FIX

    time = 0;

    for q = 1:tf
        Mz = relax_T1(M0, Mz, time, t1_vals);
        Mxy(q,:) = relax_T2s_phase_evo(sin(alpha(:,3)) .* Mz,TE_flash,t2s_vals,local_field_sim_Hz);
        Mz = cos(alpha(:,3)) .* Mz;
        time = esp;
    end
 
    M7 = Mz;
    Mxy_acq2 = Mxy;
    
    M8 = relax_T1(M0, M7, delT_M7_M8, t1_vals);
    
    % acq3
    Mz = M8;
    Mxy = zeros(tf, num_voxels, 'like', prototype_var); % FIX

    time = 0;

    for q = 1:tf
        Mz = relax_T1(M0, Mz, time, t1_vals);
        Mxy(q,:) = relax_T2s_phase_evo(sin(alpha(:,4)) .* Mz,TE_flash,t2s_vals,local_field_sim_Hz);
        Mz = cos(alpha(:,4)) .* Mz;
        time = esp;
    end
 
    M9 = Mz;
    Mxy_acq3 = Mxy;
    
    M10 = relax_T1(M0, M9, delT_M9_M10, t1_vals);
    
    % acq 4 (Multi-echo part)
    Mz = M10;
    Mxy = zeros(tf, num_voxels, 'like', prototype_var); % FIX
    Mxyecho = zeros(tf,num_voxels,length(TEs), 'like', prototype_var); % FIX: ÕâÀïµÄ10ÊÇ»Ø²¨ÊýÁ¿£¬ÐèÈ·±£·ÖÅäµÄÊÇ dlarray
    
    time = 0;

    for q = 1:tf
        Mz = relax_T1(M0, Mz, time, t1_vals);
        
        Mxy(q,:) = sin(alpha(:,5)) .* Mz;% Mxy13
        Mxy0 = sin(alpha(:,5)) .* Mz;
        
        Mxyecho(q,:,:) = relax_T2s_phase_evo(Mxy0,TEs,t2s_vals,local_field_sim_Hz);
        
        Mz = cos(alpha(:,5)) .* Mz;
        time = TR_mte;
    end
    

    Mxy_mgre0(:,:,reps) = Mxy0;
 
    M13 = Mz;

    Mstart = relax_T1(M0, M13, delT_M13_2end, t1_vals);
    

    Mz_mtx(:,:,reps) = cat(1, M1m.', M2m.', M3m.', M1.', M2.', M3.', M4.', M5.', M6.', M7.', M8.', M9.', M10.', M13.');
    
    Mxyecho_reshaped = squeeze(Mxyecho(end, :, :)).';
    Mxy_mtx(:,:,reps) = cat(1, Mxy_acq0(echo2use,:), Mxy_acq1(echo2use,:), Mxy_acq2(echo2use,:), Mxy_acq3(echo2use,:), ...
    Mxyecho_reshaped);    
    % reshape(Mxyecho(end,:,1),[1,num_voxels]),reshape(Mxyecho(end,:,2),[1,num_voxels]),...
        % reshape(Mxyecho(end,:,3),[1,num_voxels]),reshape(Mxyecho(end,:,4),[1,num_voxels]),...
        % reshape(Mxyecho(end,:,5),[1,num_voxels]),reshape(Mxyecho(end,:,6),[1,num_voxels])...
        % ,reshape(Mxyecho(end,:,7),[1,num_voxels]),reshape(Mxyecho(end,:,8),[1,num_voxels]),...
        % reshape(Mxyecho(end,:,9),[1,num_voxels]),reshape(Mxyecho(end,:,10),[1,num_voxels]));
end

end