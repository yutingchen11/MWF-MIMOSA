function write_metrics_text_loss(txtFile, run_name, cfg, feature_idx, final_train_loss, valMetrics)
% English comments only

fid = fopen(txtFile, 'w');
if fid < 0
    warning('Cannot open summary text file for writing: %s', txtFile);
    return;
end

fprintf(fid, 'Run name: %s\n', run_name);
fprintf(fid, 'Config name: %s\n', cfg.name);
fprintf(fid, 'Loss name: %s\n', cfg.loss_name);
fprintf(fid, 'Loss mode: %d\n', cfg.loss_mode);
fprintf(fid, 'Use weight: %d\n', cfg.use_weight);
fprintf(fid, 'Use norm: %d\n', cfg.use_norm);
fprintf(fid, 'epsA: %.10g\n', cfg.epsA);
fprintf(fid, 'tau: %.10g\n', cfg.tau);

fprintf(fid, '\nFeature idx: %s\n', mat2str(feature_idx));
fprintf(fid, 'Feature count: %d\n', numel(feature_idx));
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