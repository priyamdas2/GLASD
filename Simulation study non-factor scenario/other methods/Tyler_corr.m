function C_tyler = Tyler_corr(X, tol, max_iter)
%TYLER_CORR  Compute Tyler's M-estimator of correlation matrix.
%
%   C_tyler = TYLER_CORR(X) returns the Tyler correlation estimator for
%   data matrix X (n x p, rows = observations).
%
%   C_tyler = TYLER_CORR(X, tol, max_iter) uses custom tolerance and
%   maximum iteration count. Default: tol = 1e-6, max_iter = 1000.
%
%   This implementation:
%     - Uses the standard fixed-point iteration
%     - Normalizes scatter so that trace(S) = p
%     - Converts scatter to correlation (diag = 1)
%     - Applies nearestSPD at the end for numerical safety (if available)

    if nargin < 2 || isempty(tol)
        tol = 1e-6;
    end
    if nargin < 3 || isempty(max_iter)
        max_iter = 1000;
    end

    [n, p] = size(X);

    % Center columns (Tyler is scale-equivariant but not location-equivariant)
    X = X - mean(X, 1);

    % Initial scatter: sample covariance (1/n normalization)
    S = (X' * X) / n;
    % Small ridge for numerical stability
    S = S + 1e-6 * eye(p);

    for iter = 1:max_iter
        % Invert current scatter
        Sinv = inv(S);

        % Mahalanobis distances: d_i = x_i' S^{-1} x_i  (n x 1)
        XS = X * Sinv;
        d = sum(XS .* X, 2);          % each row i: x_i' S^{-1} x_i
        d = max(d, 1e-12);            % avoid division by 0

        % Weights: p / d_i
        w = p ./ d;                   % n x 1

        % Weighted scatter: (1/n) sum w_i x_i x_i'
        Xw = X .* w;                  % each row scaled by w_i
        S_new = (Xw' * X) / n;

        % Normalize trace to p (scale identification)
        S_new = S_new * (p / trace(S_new));

        % Check convergence
        rel_change = norm(S_new - S, 'fro') / max(norm(S, 'fro'), 1e-12);
        S = S_new;

        if rel_change < tol
            % fprintf('Tyler converged in %d iterations (rel change = %.2e)\n', iter, rel_change);
            break;
        end
    end

    % Convert scatter to correlation
    d_sqrt = sqrt(diag(S));
    D_inv = diag(1 ./ d_sqrt);
    C = D_inv * S * D_inv;

    % For safety, project to nearest SPD if function is available
    if exist('nearestSPD', 'file') == 2
        C = nearestSPD(C);
    end

    C_tyler = C;
end
