function [gradients, loss, loss_norm] = modelGradients_mlp_MIMOSA_lossAblation( ...
    parameters, dlXf, dlRf, alpha, loss_mode, use_weight, use_norm, epsA0, tau0)
% English comments only

if nargin < 4 || isempty(alpha)
    alpha = 0.01;
end

epsA = dlarray(single(epsA0), 'CB');
tau  = dlarray(single(tau0), 'CB');

U = mlp_model_leakyRelu(parameters, dlXf, alpha);

Ur = dlRf(1,:);  Ui = dlRf(2,:);
Pr = U(1,:);     Pi = U(2,:);

mag_gt = sqrt(Ur.^2 + Ui.^2 + epsA);

dPr = Pr - Ur;
dPi = Pi - Ui;

err_c = sqrt(dPr.^2 + dPi.^2 + epsA);

switch loss_mode

    case 1
        if use_weight
            wA = (mag_gt.^2) ./ (mag_gt.^2 + tau.^2);
        else
            wA = mag_gt .* 0 + 1;
        end

        if use_norm
            core = err_c ./ (mag_gt + epsA);
        else
            core = err_c;
        end

        loss_main = mean(wA .* core, "all");

    case 2
        loss_main = mean(0.5 * (abs(dPr) + abs(dPi)), "all");

    case 3
        loss_main = mean(0.5 * (dPr.^2 + dPi.^2), "all");

    otherwise
        error('Unknown loss_mode: %d', loss_mode);
end

loss_norm = loss_main;
loss = loss_main;

gradients = dlgradient(loss, parameters);

end