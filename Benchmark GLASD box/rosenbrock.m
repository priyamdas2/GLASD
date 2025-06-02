function y = rosenbrock(x)
% ROSENBROCK function: non-convex with narrow curved valley
% Global minimum at x = [1, 1, ..., 1], f(x) = 0

x = x(:);  % Ensure column vector
d = length(x);

% Rosenbrock sum over dimensions 1 to d-1
y = sum(100 * (x(2:d) - x(1:d-1).^2).^2 + (1 - x(1:d-1)).^2);
end