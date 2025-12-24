function [F, U] = generate_setting2_componentwise(p, r, n, nu)
% GENERATE_SETTING2_COMPONENTWISE
%   Setting 2 of Fan, Wang & Zhong (2019):
%   z_t = (f_t', u_t')' with independent coordinates, each
%   univariate t_nu, scaled so that
%       Var(f-coordinates) = 1,
%       Var(u-coordinates) = 5.
%
%   This yields Cov(z_t) = diag(I_r, 5 I_p), but the joint distribution
%   is no longer elliptical (coordinates are independent).

    d_f = r;
    d_u = p;

    if isinf(nu)
        % Gaussian limiting case: standard normals with desired variances
        % factors: N(0,1), errors: N(0,5)
        F = randn(n, d_f);                  % variance 1
        U = sqrt(5) * randn(n, d_u);        % variance 5

    else
        if nu <= 2
            warning('nu <= 2: t_nu has infinite variance; scaling to match finite variances is theoretical only.');
        end

        % Standard t_nu via normal / sqrt(chi2/nu)
        % For factors:
        g_f = randn(n, d_f);
        s_f = chi2rnd(nu, n, 1);
        t_f = g_f .* repmat(sqrt(nu ./ s_f), 1, d_f);   % n x r, standard t_nu

        % For errors:
        g_u = randn(n, d_u);
        s_u = chi2rnd(nu, n, 1);
        t_u = g_u .* repmat(sqrt(nu ./ s_u), 1, d_u);   % n x p, standard t_nu

        % Var(standard t_nu) = nu/(nu-2) for nu > 2.
        % Scale to get Var = 1 for factors and Var = 5 for u.
        if nu > 2
            baseVar = nu / (nu - 2);  % variance of t_nu
        else
            baseVar = 1;              % placeholder when variance is infinite
        end

        scale_f = sqrt(1     / baseVar);    % for Var = 1
        scale_u = sqrt(5.0   / baseVar);    % for Var = 5

        F = scale_f * t_f;      % n x r, approx Var 1 each coordinate
        U = scale_u * t_u;      % n x p, approx Var 5 each coordinate
    end
end
