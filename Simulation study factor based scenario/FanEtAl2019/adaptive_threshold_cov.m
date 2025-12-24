function S_thr = adaptive_threshold_cov(S, tau_const, n)
% ADAPTIVE_THRESHOLD_COV
%   Adaptive (entrywise) thresholding of a covariance matrix in the
%   style of Cai & Liu (2011), as used in Fan, Wang & Zhong (2019), Step 2.
%
%   For i != j:
%       tau_ij = tau * sqrt(S_ii * S_jj),
%       S_ij^T = soft_threshold(S_ij, tau_ij),
%   while diagonals are kept unchanged.
%
%   INPUT:
%       S         : p x p covariance matrix (symmetric)
%       tau_const : scalar C_tau in tau = C_tau * sqrt(log(p)/n)
%       n         : sample size
%
%   OUTPUT:
%       S_thr     : thresholded covariance matrix

    S = (S + S') / 2;
    p = size(S, 1);

    % base tau as in theory: tau ~ sqrt( log(p) / n )
    tau = tau_const * sqrt(log(p) / n);

    S_thr = zeros(p);
    diagS = diag(S);

    % keep diagonals
    S_thr(1:p+1:end) = diagS;

    % off-diagonal adaptive soft thresholding
    for i = 1:p
        for j = i+1:p
            sij = S(i, j);

            % entry-specific threshold
            d_ii = max(diagS(i), 0);
            d_jj = max(diagS(j), 0);
            tau_ij = tau * sqrt(d_ii * d_jj);

            val = soft_threshold(sij, tau_ij);

            S_thr(i, j) = val;
            S_thr(j, i) = val;
        end
    end
end

function y = soft_threshold(x, t)
% SOFT_THRESHOLD
%   y = sign(x) * max(|x| - t, 0)

    if t <= 0
        y = x;
        return;
    end

    y = sign(x) .* max(abs(x) - t, 0);
end
