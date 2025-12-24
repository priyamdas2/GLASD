function mu = huber_mean_newton(w, alpha, max_iter, tol)
% HUBER_MEAN_NEWTON
%   1D Huber M-estimator for the mean using Newton iterations.
%
%   Solve for mu in:
%       sum_t psi_alpha(w_t - mu) = 0,
%   where psi_alpha is the Huber score function:
%       psi(u) = u                   if |u| <= alpha
%              = alpha * sign(u)     if |u| >  alpha
%
%   INPUT:
%       w        : n x 1 data vector
%       alpha    : Huber threshold (>0)
%       max_iter : maximum number of Newton iterations
%       tol      : stopping tolerance on parameter change
%
%   OUTPUT:
%       mu       : robust mean estimate

    if nargin < 3, max_iter = 50;  end
    if nargin < 4, tol      = 1e-6; end

    w = w(:);
    n = numel(w);

    % initialization at sample mean
    mu = mean(w);

    for it = 1:max_iter
        u = w - mu;

        % psi_alpha
        psi = u;
        idx_big = abs(u) > alpha;
        psi(idx_big) = alpha .* sign(u(idx_big));

        g = sum(psi);  % gradient: should be 0 at optimum

        % derivative wrt mu:
        % d/dmu psi(w_t - mu) = -1 for |w_t - mu| <= alpha, 0 otherwise
        idx_in = abs(u) <= alpha;
        gprime = -sum(idx_in);

        if gprime == 0
            % no inliers under current mu; break or take tiny step
            break;
        end

        step = g / gprime;        % Newton step
        mu_new = mu + step;

        if abs(mu_new - mu) < tol * max(1, abs(mu))
            mu = mu_new;
            return;
        end

        mu = mu_new;
    end
end
