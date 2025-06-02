function [C_best, f_best, history] = GLASD_PD(f, C0, params)
% GLASD: Global Adaptive Stochastic Descent on PD matrix with unit
% diagonals
%
% Inputs:
%   f      - function handle to minimize (on PD matrix with unit diagonals)
%   C0     - initial point, PD matrix with unit diagonals
%   params - (optional) struct of algorithm parameters
%
% Outputs:
%   C_best  - best point in original domain
%   f_best  - best function value found
%   history - struct with .fvals and .C_best (best in scaled space)

if ~isequal(C0, C0') || any(eig(C0) < 0)
    error('C0 must be a symmetric positive semi-definite matrix');
end
if any(abs(diag(C0) - 1) > 1e-10)
    error('C0 must have 1s on its diagonal (within tolerance)');
end

x0 = Corr2Theta(C0);
n = length(x0);
[lb, ub] = angular_bounds(size(C0,1));
lb = lb(:); ub = ub(:); x0 = x0(:);
range = ub - lb;
z0 = (x0 - lb) ./ range;

if any(z0 < 0 | z0 > 1)
    error('All elements of z0 must lie within [0, 1]');
end

% Scaled objective function
f_scaled_of_C = @(C) f(C);

% --- Default parameters ---
default.s_init  = 0.1;
default.s_inc   = 2;
default.s_dec   = 2;
default.p_inc   = 2;
default.p_dec   = 2;
default.m       = 5;
default.c       = 0.001*log(n);
default.T       = round(3000 * log(n));
default.M       = 4*n;
default.epsilon = 1e-20;

if nargin < 5
    params = default;
else
    fnames = fieldnames(default);
    for k = 1:length(fnames)
        if ~isfield(params, fnames{k})
            params.(fnames{k}) = default.(fnames{k});
        end
    end
end

% --- Initialization in scaled space ---
z = z0;
p_init  = (1/(2*n))*ones(2*n,1);
s = params.s_init * ones(2*n,1);
p = p_init(:); p = p / sum(p);
C_z = Theta2Corr(lb + range .* z);
f_curr = f_scaled_of_C(C_z);
f_best = f_curr;
z_best = z;
fvals = f_curr;
f_window = f_best * ones(params.M, 1);

for t = 1:params.T
    if rand < 1 - 1/params.m
        % Greedy descent
        j = randsample(2*n, 1, true, p);
        i = ceil(j/2);
        sign = 1 - 2 * mod(j,2);  % +1 if odd, -1 if even
        if sign > 0
            delta_i = min(s(j), (1 - z(i)) / 2);
        else
            delta_i = -min(s(j), z(i) / 2);
        end
        delta = zeros(n,1); delta(i) = delta_i;
        z_new = z + delta;
        C_z_new = Theta2Corr(lb + range .* z_new);
        f_new = f_scaled_of_C(C_z_new);
        if f_new < f_curr
            z = z_new; f_curr = f_new;
            s(j) = s(j) * params.s_inc;
            p(j) = p(j) * params.p_inc;
            p = p / sum(p);
        else
            s(j) = s(j) / params.s_dec;
            p(j) = p(j) / params.p_dec;
            p = p / sum(p);
        end
    else
        % Forced exploration
        i = randi(n);
        sign = 2 * randi(2) - 3;
        if sign > 0
            delta_i = rand * (1 - z(i));
        else
            delta_i = -rand * z(i);
        end
        delta = zeros(n,1); delta(i) = delta_i;
        z_new = z + delta;
        C_z_new = Theta2Corr(lb + range .* z_new);
        f_new = f_scaled_of_C(C_z_new);
        q_t = min(1, params.m * params.c / log(1 + t));
        if f_new < f_curr || rand < q_t
            z = z_new; f_curr = f_new;
        end
    end

    % Update best solution
    if f_curr < f_best
        f_best = f_curr;
        z_best = z;
        C_best = Theta2Corr(z_best);
    end
    f_window = [f_window(2:end); f_best];
    fvals(end+1) = f_best;

    if t >= params.M && f_window(1) - f_best < params.epsilon
        break;
    end
end

x_best = lb + range .* z_best;  % Transform back to original space
C_best = Theta2Corr(x_best);
history.fvals = fvals;
history.C_best = C_best;
end
