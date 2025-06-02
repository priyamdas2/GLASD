function y = sumsquares(x)
% SUMSQUARES function: simple convex function with global minimum at 0
% Global minimum at x = [0, 0, ..., 0], f(x) = 0

x = x(:);  % Ensure x is a column vector
d = length(x);
coeff = (1:d)';  % Weights: 1, 2, ..., d

y = sum(coeff .* (x .^ 2));
end