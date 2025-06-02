function y = ackley(x)
% ACKLEY function: multimodal test function
% Global minimum at x = 0, f(x) = 0
a = 20;
b = 0.2;
c = 2*pi;
d = length(x);

sum1 = sum(x.^2);
sum2 = sum(cos(c * x));
term1 = -a * exp(-b * sqrt(sum1 / d));
term2 = -exp(sum2 / d);
y = term1 + term2 + a + exp(1);
end