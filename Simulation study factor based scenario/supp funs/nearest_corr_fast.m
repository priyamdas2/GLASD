function C0 = nearest_corr_fast(C, eps)

    if nargin < 2
        eps = 1e-6;
    end

    % 1. Symmetrize
    C = (C + C')/2;

    % 2. Convert to correlation scale
    D = diag(1./sqrt(diag(C)));
    R = D*C*D;

    % 3. Eigenvalue clipping
    [V, D_eig] = eig(R);
    d = max(diag(D_eig), eps);
    R_pd = V*diag(d)*V';

    % 4. Renormalize to enforce exact unit diagonal
    D2 = diag(1./sqrt(diag(R_pd)));
    C0 = D2*R_pd*D2;

    % 5. Final numerical cleanup
    C0 = (C0 + C0')/2;
end
