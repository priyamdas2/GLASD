function y = rastrigin(x)
% RASTRIGIN function: highly multimodal with large search space
% Global minimum at x = 0, f(x) = 0

A = 10;
d = length(x);
y = A * d + sum(x.^2 - A * cos(2 * pi * x));
end