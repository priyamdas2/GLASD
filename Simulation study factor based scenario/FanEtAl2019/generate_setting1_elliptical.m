function [F, U] = generate_setting1_elliptical(p, r, n, nu)
% GENERATE_SETTING1_ELLIPTICAL
%   Setting 1 of Fan, Wang & Zhong (2019):
%   z_t = (f_t', u_t')' from multivariate t_nu with
%   Cov(z_t) = diag(I_r, 5 I_p).
%
%   We implement this via the standard multivariate t construction:
%       g ~ N(0, Sigma0), s ~ chi2_nu, independent,
%       z = g ./ sqrt(s / nu),
%   where Cov(z) = (nu/(nu-2)) * Sigma0 for nu > 2.
%   To get Cov(z) = diag(I_r, 5I_p), choose:
%       Sigma0 = ((nu-2)/nu) * diag(I_r, 5I_p).

    d = r + p;  % total dimension of z_t

    if isinf(nu)
        % Gaussian case: z_t ~ N(0, diag(I_r, 5 I_p))
        scale_diag = [ones(r, 1); 5 * ones(p, 1)];     % diag entries
        std_vec    = sqrt(scale_diag);                 % std deviations

        Z = randn(n, d) .* repmat(std_vec', n, 1);     % n x d

    else
        if nu <= 2
            warning('nu <= 2: theoretical covariance of t_nu is infinite; using same construction anyway.');
        end

        % Desired covariance of z_t
        Cov_z = blkdiag(eye(r), 5 * eye(p));           % diag(I_r, 5 I_p)

        % Under multivariate t_nu, Cov(z) = (nu/(nu-2)) * Sigma0.
        % So pick Sigma0 = ((nu-2)/nu) * Cov_z.
        factor = (nu - 2) / nu;
        Sigma0 = factor * Cov_z;

        % Cholesky factor of Sigma0
        % (since Sigma0 is diagonal, this is simple, but we keep it general)
        A = chol(Sigma0, 'lower');                     % d x d

        % draw g ~ N(0, Sigma0)
        G = randn(n, d) * A';                          % n x d

        % chi-square draws
        s = chi2rnd(nu, n, 1);                         % n x 1

        % multivariate t draws: each row scaled by sqrt(nu./s)
        scale = sqrt(nu ./ s);                         % n x 1
        Z = G .* repmat(scale, 1, d);                  % n x d
    end

    % split into factors and idiosyncratic components
    F = Z(:, 1:r);         % n x r
    U = Z(:, r+1:end);     % n x p
end
