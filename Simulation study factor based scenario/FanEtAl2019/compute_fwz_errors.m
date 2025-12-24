%% ============================================
% Helper: function to fit FWZ and compute errors
% ============================================

function err = compute_fwz_errors(Y, F, Sigma_true, Sigma_u_true, Sigma_inv_true)
    % Fit Fan–Wang–Zhong robust covariance estimator
    opts = struct();
    opts.tau_const = 2;   % tau = 2 * sqrt(log(p)/n) as in paper
    % other options (alpha multipliers, etc.) use defaults inside fan2019_cov_estimator

    fit = fan2019_cov_estimator(Y, F, opts);

    Sigma_hat   = fit.Sigma_hat;
    Sigma_u_hat = fit.Sigma_u_hat;

    % Basic sanity symmetrization (should already be symmetric)
    Sigma_hat   = (Sigma_hat   + Sigma_hat') / 2;
    Sigma_u_hat = (Sigma_u_hat + Sigma_u_hat') / 2;

    % Ensure invertible (FWZ should give PD, but just in case)
    [Vh, Dh] = eig((Sigma_hat + Sigma_hat') / 2);
    lam      = diag(Dh);
    eps_reg  = 1e-6 * max(lam);
    Dh_reg   = diag(max(lam, eps_reg));
    Sigma_hat_reg = Vh * Dh_reg * Vh';

    Sigma_inv_hat = inv(Sigma_hat_reg);

    % Frobenius norm errors
    err_Sigma_frob    = norm(Sigma_hat   - Sigma_true,   'fro');
    err_Sigma_op      = norm(Sigma_hat   - Sigma_true);        % operator (2-)norm
    err_Prec_frob     = norm(Sigma_inv_hat - Sigma_inv_true, 'fro');
    err_Prec_op       = norm(Sigma_inv_hat - Sigma_inv_true); % operator norm
    err_Sigma_u_frob  = norm(Sigma_u_hat - Sigma_u_true, 'fro');
    err_Sigma_u_op    = norm(Sigma_u_hat - Sigma_u_true);

    % Package as a row vector
    err = [err_Sigma_frob, err_Sigma_op, ...
           err_Prec_frob,  err_Prec_op, ...
           err_Sigma_u_frob, err_Sigma_u_op];
end
