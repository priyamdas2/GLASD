function val = loss_tukey(X, C)
    % LOSS_TUKEY computes Tukey biweight-based robust loss
    % Input:
    %   X : n x p data matrix
    %   C : p x p positive definite covariance (or correlation) matrix
    % Output:
    %   val : scalar loss value
    
    [n, p] = size(X);
    C_inv = inv(C);
    d = mahal_squared_distances(X, C_inv);  % Mahalanobis squared distances

    % Tukey biweight threshold (chi-squared based)
    threshold = quantile(d, 0.75) + 3 * iqr(d);
    c = sqrt(threshold);  % scale parameter (applied to sqrt Mahalanobis dist)

    % Apply Tukey biweight loss to sqrt distances
    r = sqrt(d);  % Mahalanobis distances
    rho = zeros(n, 1);
    
    inside = r <= c;
    r_scaled = r(inside) / c;
    rho(inside) = (c^2 / 6) * (1 - (1 - r_scaled.^2).^3);
    rho(~inside) = c^2 / 6;

    % Final loss: log(det(C)) + sum of robustified distances
    val = 0.5 * n * log(det(C)) + sum(rho);
end
