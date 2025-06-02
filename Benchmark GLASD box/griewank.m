function y = griewank(x)
% GRIEWANK function: multimodal with many regularly spaced local minima
% Global minimum at x = 0, f(x) = 0

d = length(x);
sum_term = sum(x.^2) / 4000;
prod_term = prod(cos(x ./ sqrt((1:d)')));
y = sum_term - prod_term + 1;
end