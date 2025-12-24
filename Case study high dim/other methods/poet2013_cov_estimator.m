function out = poet2013_cov_estimator(Y, r, opts)
% POET2013_COV_ESTIMATOR
% Classical POET estimator from:
% Fan, Liao & Mincheva (2013, JRSS-B)
% Large covariance estimation by thresholding principal orthogonal complements
%
% Model: y_t = B f_t + u_t, with factors UNOBSERVED
%
% INPUT:
%   Y    : n x p data matrix
%   r    : number of factors (K)
%   opts : options structure (optional)
%
% OUTPUT (same naming as FWZ):
%   out.Sigma_hat
%   out.Sigma_u_hat
%   out.Sigma_u_raw
%   out.Sigma11_hat   (sample covariance)
%   out.Sigma12_hat   (factor part)
%   out.Sigma22_hat   (factor covariance)
%   out.F_hat         (estimated factors)
%   out.B_hat         (estimated loadings)

    if nargin < 3
        opts = struct();
    end

    % ---- defaults ----
    if ~isfield(opts, 'center'),   opts.center   = true; end
    if ~isfield(opts, 'tau_const'),opts.tau_const= 2.0;  end
    if ~isfield(opts, 'threshold'),opts.threshold='hard'; end

    [n, p] = size(Y);

    % ---- Step 0: center data ----
    if opts.center
        Y = Y - mean(Y, 1);
    end

    % ---- Step 1: sample covariance ----
    Sigma_sample = (Y' * Y) / n;
    Sigma_sample = (Sigma_sample + Sigma_sample') / 2;

    % ---- Step 2: PCA for latent factors ----
    [V, D] = eig(Sigma_sample);
    [lambda, idx] = sort(diag(D), 'descend');
    V = V(:, idx);

    V_r = V(:, 1:r);
    Lambda_r = diag(lambda(1:r));

    % Estimated factor covariance part
    Sigma_lowrank = V_r * Lambda_r * V_r';

    % ---- Step 3: Raw idiosyncratic covariance ----
    Sigma_u_raw = Sigma_sample - Sigma_lowrank;
    Sigma_u_raw = (Sigma_u_raw + Sigma_u_raw') / 2;

    % ---- Step 4: Thresholding Sigma_u ----
    tau = opts.tau_const * sqrt(log(p) / n);
    Sigma_u_hat = Sigma_u_raw;

    for i = 1:p
        for j = 1:p
            if i ~= j
                if strcmp(opts.threshold, 'hard')
                    if abs(Sigma_u_hat(i,j)) < tau
                        Sigma_u_hat(i,j) = 0;
                    end
                elseif strcmp(opts.threshold, 'soft')
                    Sigma_u_hat(i,j) = sign(Sigma_u_hat(i,j)) * max(abs(Sigma_u_hat(i,j)) - tau, 0);
                end
            end
        end
    end
    Sigma_u_hat = (Sigma_u_hat + Sigma_u_hat') / 2;

    % ---- Step 5: Final covariance ----
    Sigma_hat = Sigma_lowrank + Sigma_u_hat;
    Sigma_hat = (Sigma_hat + Sigma_hat') / 2;

    % ---- Recover factors and loadings (optional but useful) ----
    F_hat = Y * V_r;              % n x r
    B_hat = V_r * sqrt(Lambda_r);

    % ---- Package outputs (FWZ-compatible) ----
    out = struct();
    out.Sigma_hat   = Sigma_hat;
    out.Sigma_u_hat = Sigma_u_hat;
    out.Sigma_u_raw = Sigma_u_raw;
    out.Sigma11_hat = Sigma_sample;
    out.Sigma12_hat = Sigma_lowrank;
    out.Sigma22_hat = Lambda_r;
    out.F_hat       = F_hat;
    out.B_hat       = B_hat;
end
