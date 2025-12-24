function Ahat = nearestSPD(A)
% NEARESTSPD - Find the nearest Symmetric Positive Definite matrix to A
%   Ahat = nearestSPD(A)
%
%   This function returns the nearest Symmetric Positive Definite matrix
%   to a given input matrix A using the method of Higham (1988).

    % Make symmetric
    B = (A + A') / 2;

    % Eigenvalue decomposition
    [V, D] = eig(B);
    D = diag(D);

    % Replace negative eigenvalues with small positive values
    D(D < 0) = 0;
    Ahat = V * diag(D) * V';

    % Ensure symmetry again
    Ahat = (Ahat + Ahat') / 2;

    % Add a small jitter if needed to ensure positive definiteness
    k = 0;
    while ~isPositiveDefinite(Ahat)
        k = k + 1;
        min_eig = min(eig(Ahat));
        Ahat = Ahat + eye(size(A)) * (-min_eig * (1 + 1e-6));
        if k > 5
            warning('nearestSPD: more than 5 adjustments made to ensure SPD.');
            break;
        end
    end
end

function tf = isPositiveDefinite(M)
    [~, p] = chol(M);
    tf = (p == 0);
end
