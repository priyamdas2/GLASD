function out = fan2019_cov_estimator_unobservedF(Y, r, opts)
% FAN2019_COV_ESTIMATOR_UNOBSERVEDF
% Robust covariance estimation for approximate factor models
% following Fan, Wang & Zhong (2019, J. Econometrics).
%
% Model: y_t = B f_t + u_t, with factors F UNOBSERVED.
% Factors are estimated via robust PCA as in FWZ Section 4.
%
% INPUT:
%   Y    : n x p data matrix
%   r    : number of factors
%
% OUTPUT (same fields as observed-F version):
%   .Sigma_hat
%   .Sigma_u_hat
%   .Sigma_u_raw
%   .Sigma_z_hat   (here equals Sigma_Y_hat)
%   .Sigma11_hat
%   .Sigma12_hat
%   .Sigma22_hat

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

    [n, p] = size(Y);

    % ---- Step 0: center data ----
    if opts.center
        Y = Y - mean(Y, 1);
    end

    %% ============================================================
    % Step 1: Robust covariance of Y only (Huber entrywise)
    % ============================================================

    Sigma_Y_hat = zeros(p, p);
    log_term = log(p);

    for i = 1:p
        yi = Y(:, i);
        for j = i:p
            yj = Y(:, j);
            w  = yi .* yj;

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
            Sigma_Y_hat(i, j) = mu_hat;
            Sigma_Y_hat(j, i) = mu_hat;
        end
    end
    Sigma_Y_hat = (Sigma_Y_hat + Sigma_Y_hat') / 2;

    %% ============================================================
    % Step 2: Robust PCA to estimate factor space
    % ============================================================

    [V, D] = eig(Sigma_Y_hat);
    [lam, idx] = sort(diag(D), 'descend');
    V = V(:, idx);

    % Estimated factor loadings
    B_hat = V(:, 1:r) * diag(sqrt(lam(1:r)));

    % Low-rank component
    Sigma_hat_lowrank = B_hat * B_hat';

    %% ============================================================
    % Step 3: Raw idiosyncratic covariance
    % ============================================================

    Sigma_u_raw = Sigma_Y_hat - Sigma_hat_lowrank;
    Sigma_u_raw = (Sigma_u_raw + Sigma_u_raw') / 2;

    %% ============================================================
    % Step 4: Adaptive thresholding (FWZ universal τ)
    % ============================================================

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

    %% ============================================================
    % Step 5: Final covariance estimator
    % ============================================================

    Sigma_hat = Sigma_hat_lowrank + Sigma_u_hat;
    Sigma_hat = (Sigma_hat + Sigma_hat') / 2;

    %% ============================================================
    % Package outputs (SAME field names as observed-F version)
    % ============================================================

    out = struct();
    out.Sigma_hat   = Sigma_hat;
    out.Sigma_u_hat = Sigma_u_hat;
    out.Sigma_u_raw = Sigma_u_raw;
    out.Sigma_z_hat = Sigma_Y_hat;       % now just Σ_Y
    out.Sigma11_hat = Sigma_Y_hat;       % for compatibility
    out.Sigma12_hat = B_hat;
    out.Sigma22_hat = eye(r);
end
