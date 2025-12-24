function out = fwz_generate_data(p, r, n, nu, settingType, B)
% FWZ_GENERATE_DATA
%   Data generation scheme mimicking Fan, Wang & Zhong (2019, Sec. 4.2).
%
%   We generate n samples of z_t = (f_t', u_t')' under two settings:
%
%   Setting 1 ('elliptical'):
%       z_t ~ multivariate t_nu(0, Cov_z), where
%       Cov_z = diag(I_r, 5 * I_p).
%       This is an elliptical multivariate t with the stated covariance.
%
%   Setting 2 ('componentwise'):
%       Each coordinate of z_t is independent univariate t_nu, scaled so that
%       Var(f-coordinates) = 1, Var(u-coordinates) = 5, hence
%       Cov_z = diag(I_r, 5 * I_p), but the joint distribution is non-elliptical.
%
%   In both settings, we generate a loading matrix B (p x r), with
%   each row B(i,:) ~ N_r(0, I_r). Once generated, B is treated as fixed.
%
%   Model:
%       y_t = B f_t + u_t,   t = 1,...,n
%
%   INPUT:
%       p           : dimension of y_t (number of observed variables)
%       r           : number of factors
%       n           : sample size
%       nu          : degrees of freedom (e.g. 3, 5, or Inf for Gaussian)
%       settingType : 'elliptical' or 'componentwise'
%       B (optional): p x r loading matrix. If omitted or empty, generated here.
%
%   OUTPUT (struct):
%       .Y         : n x p matrix of responses (rows = y_t')
%       .F         : n x r matrix of factors (rows = f_t')
%       .U         : n x p matrix of idiosyncratic errors (rows = u_t')
%       .B         : p x r loading matrix used
%       .Sigma_u   : p x p true Σ_u = 5 * I_p
%       .Sigma_y   : p x p true Σ = B B' + 5 * I_p
%       .Cov_z     : (r+p) x (r+p) covariance of z_t (should be diag(I_r, 5I_p))

    % ---- basic checks ----
    if nargin < 5 || isempty(settingType)
        error('You must specify settingType = ''elliptical'' or ''componentwise''.');
    end

    if nargin < 6
        B = [];
    end

    % ---- generate B if not supplied ----
    if isempty(B)
        % Each row of B ~ N_r(0, I_r)
        B = randn(p, r);
    else
        % light check
        if ~isequal(size(B), [p, r])
            error('Provided B must be of size p x r.');
        end
    end

    % true idiosyncratic covariance and total covariance
    Sigma_u = 5 * eye(p);
    Sigma_y = B * B' + Sigma_u;

    % ---- generate (F, U) under the chosen setting ----
    switch lower(settingType)
        case 'elliptical'
            [F, U] = generate_setting1_elliptical(p, r, n, nu);

        case 'componentwise'
            [F, U] = generate_setting2_componentwise(p, r, n, nu);

        otherwise
            error('Unknown settingType. Use ''elliptical'' or ''componentwise''.');
    end

    % responses
    Y = F * B' + U;   % n x p

    % pack outputs
    out = struct();
    out.Y       = Y;
    out.F       = F;
    out.U       = U;
    out.B       = B;
    out.Sigma_u = Sigma_u;
    out.Sigma_y = Sigma_y;
    out.Cov_z   = blkdiag(eye(r), 5 * eye(p));  % target Cov(z_t)
end
