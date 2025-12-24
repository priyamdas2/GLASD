function [x_best, f_best, history] = GLASD(f, lb, ub, x0, params)
% GLASD: Global Adaptive Stochastic Descent with [0,1]^n scaling
%
% Inputs:
%   f      - function handle to minimize (on original scale)
%   lb     - lower bounds (n x 1)
%   ub     - upper bounds (n x 1)
%   x0     - initial point (n x 1)
%   params - (optional) struct of algorithm parameters
%
% Outputs:
%   x_best  - best point in original domain
%   f_best  - best function value found
%   history - struct with .fvals and .z_best (best in scaled space)

n = length(x0);
if length(lb) ~= n || length(ub) ~= n
    error('lb, ub, and x0 must have the same length.');
end

lb = lb(:); ub = ub(:); x0 = x0(:);
range = ub - lb;
z0 = (x0 - lb) ./ range;

% Scaled objective function
f_scaled = @(z) f(lb + range .* z);

% --- Default parameters ---
default.s_init  = 0.5;
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
f_curr = f_scaled(z);
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
        f_new = f_scaled(z_new);
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
        f_new = f_scaled(z_new);
        q_t = min(1, params.m * params.c / log(1 + t));
        if f_new < f_curr || rand < q_t
            z = z_new; f_curr = f_new;
        end
    end

    % Update best solution
    if f_curr < f_best
        f_best = f_curr;
        z_best = z;
    end
    f_window = [f_window(2:end); f_best];
    fvals(end+1) = f_best;

    if t >= params.M && f_window(1) - f_best < params.epsilon
        break;
    end
end

x_best = lb + range .* z_best;  % Transform back to original space
history.fvals = fvals;
history.z_best = z_best;
end
