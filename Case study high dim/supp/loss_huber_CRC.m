function val = loss_huber_CRC(X, C)
    [n, p] = size(X);
    C_inv = inv(C);
    d = mahal_squared_distances(X, C_inv);
    d_sqrt = sqrt(d);
    
    delta = quantile(d, 0.90) + 3 * iqr(d);
    sqrt_delta = sqrt(delta);
    
    h = zeros(size(d));
    mask = d <= delta;
    h(mask) = d(mask);            % Huber for squared distances
    h(~mask) = 2*sqrt_delta*d_sqrt(~mask) - delta;  % linear region

    val = 0.5*sum(h) + 0.5 * n * log(det(C));
end