function out = fan2019_cov_estimator(Y, F, opts)
% FAN2019_COV_ESTIMATOR
% Robust covariance estimation for approximate factor models
% following Fan, Wang & Zhong (2019, J. Econometrics).
%
% Model: y_t = B f_t + u_t, with factors F observed (known).
%
% OUTPUT:
%   Same fields as your original version (unchanged)

    if nargin < 3
        opts = struct();
    end

    % ---- defaults ----
    if ~isfield(opts, 'center'),          opts.center          = true;   end
    if ~isfield(opts, 'alpha_mult_diag'), opts.alpha_mult_diag = 1.0;    end
    if ~isfield(opts, 'alpha_mult_off'),  opts.alpha_mult_off  = 0.8;    end
    if ~isfield(opts, 'tau_const'),       opts.tau_const       = 2.0;    end
    if ~isfield(opts, 'max_iter'),        opts.max_iter        = 5000;     end
    if ~isfield(opts, 'tol'),             opts.tol             = 1e-6;   end

    [n, p]  = size(Y);
    [n2, r] = size(F);
    if n2 ~= n
        error('Y and F must have the same number of rows.');
    end

    % ---- Step 0: center data ----
    if opts.center
        Y = Y - mean(Y, 1);
        F = F - mean(F, 1);
    end

    % Enforce Cov(F) = I_r (FWZ assumption)
    F = zscore(F);

    % concatenate z_t = [y_t, f_t]
    Z = [Y, F];            % n x (p+r)
    d = p + r;

    % ---- Step 1: entrywise robust covariance of Z ----
    Sigma_z_hat = zeros(d, d);
    log_term = log(d);

    for i = 1:d
        zi = Z(:, i);
        for j = i:d
            zj = Z(:, j);
            w  = zi .* zj;

            % robust scale via MAD
            med_w = median(w);
            rob_scale = median(abs(w - med_w)) / 0.6745;
            if rob_scale <= 0 || ~isfinite(rob_scale)
                rob_scale = std(w);
            end
            if rob_scale <= 0 || ~isfinite(rob_scale)
                rob_scale = 1;
            end

            % FWZ-style alpha
            if i == j
                alpha = opts.alpha_mult_diag * rob_scale * sqrt(n / log_term);
            else
                alpha = opts.alpha_mult_off  * rob_scale * sqrt(n / log_term);
            end

            mu_hat = huber_mean_newton(w, alpha, opts.max_iter, opts.tol);
            Sigma_z_hat(i, j) = mu_hat;
            Sigma_z_hat(j, i) = mu_hat;
        end
    end

    % extract blocks
    Sigma11_hat = Sigma_z_hat(1:p,      1:p);
    Sigma12_hat = Sigma_z_hat(1:p,      p+1:end);
    Sigma21_hat = Sigma_z_hat(p+1:end,  1:p);
    Sigma22_hat = Sigma_z_hat(p+1:end,  p+1:end);

    Sigma22_hat = (Sigma22_hat + Sigma22_hat') / 2;
    Sigma22_inv = inv(Sigma22_hat);

    % ---- Step 2: raw Sigma_u ----
    Sigma_u_raw = Sigma11_hat - Sigma12_hat * Sigma22_inv * Sigma21_hat;
    Sigma_u_raw = (Sigma_u_raw + Sigma_u_raw') / 2;

    % ---- Step 2b: universal threshold τ = 2*sqrt(log p/n) ----
    tau = opts.tau_const * sqrt(log(p) / n);
    Sigma_u_hat = Sigma_u_raw;
    for i = 1:p
        for j = 1:p
            if i ~= j && abs(Sigma_u_hat(i,j)) < tau
                Sigma_u_hat(i,j) = 0;
            end
        end
    end
    Sigma_u_hat = (Sigma_u_hat + Sigma_u_hat') / 2;

    % ---- Step 3: final covariance reconstruction ----
    Sigma_hat_lowrank = Sigma12_hat * Sigma22_inv * Sigma21_hat;
    Sigma_hat = Sigma_hat_lowrank + Sigma_u_hat;
    Sigma_hat = (Sigma_hat + Sigma_hat') / 2;

    % ---- package outputs (unchanged names) ----
    out = struct();
    out.Sigma_hat   = Sigma_hat;
    out.Sigma_u_hat = Sigma_u_hat;
    out.Sigma_u_raw = Sigma_u_raw;
    out.Sigma_z_hat = Sigma_z_hat;
    out.Sigma11_hat = Sigma11_hat;
    out.Sigma12_hat = Sigma12_hat;
    out.Sigma22_hat = Sigma22_hat;
end
