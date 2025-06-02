function val = loss_truncated(X, C)
    [n, p] = size(X);
    C_inv = inv(C);
    d = mahal_squared_distances(X, C_inv);
    
%     Q1 = prctile(d, 25);
%     Q3 = prctile(d, 75);
%     IQR = Q3 - Q1;
%     threshold = Q3 + 1.5 * IQR;
%     
    threshold = quantile(d, 0.75) + 3 * iqr(d);
    outlier_flags = d > threshold;

    % Truncate: keep distances for non-outliers, cap outliers at max non-outlier
    d_cap = d;
    d_cap(outlier_flags) = threshold;

    % Compute objective
    val = 0.5 * n * log(det(C)) + 0.5 * sum(d_cap);
end
