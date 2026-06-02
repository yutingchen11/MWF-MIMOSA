
function seq = build_MIMOSA_acq_times(seq, echo2use)
% English comments only
% Build acquisition-time features aligned with Mxy_mtx 14 points.
% Supports echo2use and acq0/acq1 initial time offset (9.7ms) used in simulator.
%
% Inputs:
%   seq      : struct with required fields
%   echo2use : scalar echo index used for acq0..acq3 (1-based)
%
% Outputs:
%   seq.t_acq14_s      : 1x14 absolute acquisition times within one repetition (seconds)
%   seq.t_sinceInv14_s : 1x14 time since inversion (seconds), zero for pre-inversion points

% ---- Defaults ----
if nargin < 2 || isempty(echo2use)
    echo2use = 1;
end
if ~isscalar(echo2use) || echo2use < 1 || floor(echo2use) ~= echo2use
    error('build_MIMOSA14_acq_times:BadEchoIndex', 'echo2use must be a positive integer.');
end

% ---- Validate required fields ----
reqFields = {'esp','turbo_factor','TR_mte','TEs','gap_between_readouts','TE_flash'};
for k = 1:numel(reqFields)
    if ~isfield(seq, reqFields{k})
        error('build_MIMOSA14_acq_times:MissingField', ...
            'seq.%s is required.', reqFields{k});
    end
end

esp = double(seq.esp);
tf  = double(seq.turbo_factor);
etl = tf * esp;

TE_flash = double(seq.TE_flash);
TEs_mte  = double(seq.TEs(:).');     % 1x10
TR_mte   = double(seq.TR_mte);

necho = length(TEs_mte);

% ---- Simulator timing constants ----
delT_M1m_M2m = 29.7e-3;
delT_M0_M1m  = double(seq.gap_between_readouts) - etl - delT_M1m_M2m;

delT_M3m_M1  = delT_M0_M1m;
delT_M1_M2   = 89.7e-3;

delT_M5_M6   = 100e-3 - 6.45e-3;
delT_M7_M8   = 1e-3;
delT_M9_M10  = 1e-3;

delT_M2_M6 = double(seq.gap_between_readouts);
delT_M2_M3 = etl;
delT_M4_M5 = 12.8e-3;
delT_M3_M4 = delT_M2_M6 - delT_M2_M3 - delT_M4_M5 - delT_M5_M6;

% ---- Echo time offsets within FLASH readout ----
% In sim, acq0/acq1 use time=9.7ms before first echo sampling; acq2/acq3 use time=0.
time0_acq01 = 9.7e-3;
time0_acq23 = 0;

% Sampling offset for selected echo within a FLASH readout:
% relax_T1(..., time0) happens BEFORE sampling, then TE_flash is applied in relax_T2s_phase_evo.
t_echo_offset_acq01 = time0_acq01 + (double(echo2use) - 1) * esp + TE_flash;
t_echo_offset_acq23 = time0_acq23 + (double(echo2use) - 1) * esp + TE_flash;

% End-of-readout duration approximation (for chaining blocks)
% In sim, each q uses time=esp after the first, but we approximate block span with etl.
% This is sufficient for feature timing consistency across points.
readout_span = etl;

% ---- Build absolute acquisition times for each of 14 Mxy points ----
t_acq14 = zeros(1,4+necho);

% acq0 starts after delT_M0_M1m + delT_M1m_M2m
t_start_acq0 = delT_M0_M1m + delT_M1m_M2m;
t_acq14(1)   = t_start_acq0 + t_echo_offset_acq01;
t_end_acq0   = t_start_acq0 + readout_span;

% acq1 starts after end of acq0 + delT_M3m_M1 + delT_M1_M2
t_start_acq1 = t_end_acq0 + delT_M3m_M1 + delT_M1_M2;
t_acq14(2)   = t_start_acq1 + t_echo_offset_acq01;
t_end_acq1   = t_start_acq1 + readout_span;

% inversion happens at M5 right after reaching M4
t_inv = t_end_acq1 + delT_M3_M4;

% acq2 starts after inversion + delT_M5_M6
t_start_acq2 = t_inv + delT_M5_M6;
t_acq14(3)   = t_start_acq2 + t_echo_offset_acq23;
t_end_acq2   = t_start_acq2 + readout_span;

% acq3 starts after end of acq2 + delT_M7_M8
t_start_acq3 = t_end_acq2 + delT_M7_M8;
t_acq14(4)   = t_start_acq3 + t_echo_offset_acq23;
t_end_acq3   = t_start_acq3 + readout_span;

% multi-echo starts after end of acq3 + delT_M9_M10
t_start_mte = t_end_acq3 + delT_M9_M10;

% store Mxyecho(end,:,e): q=tf => (tf-1)*TR_mte + TE(e)
if numel(TEs_mte) ~= necho
    error('build_MIMOSA14_acq_times:TECount', ...
        'Expected seq.TEs to have 10 echoes for the multi-echo block, got %d.', numel(TEs_mte));
end
t_qend = (tf - 1) * TR_mte;
for e = 1:necho
    t_acq14(4+e) = t_start_mte + t_qend + TEs_mte(e);
end

% ---- Time since inversion (0 for pre-inversion points) ----
t_sinceInv14 = max(t_acq14 - t_inv, 0);

% ---- Cast to single for minibatch consistency ----
seq.t_acq14_s      = single(t_acq14);
seq.t_sinceInv14_s = single(t_sinceInv14);

end
