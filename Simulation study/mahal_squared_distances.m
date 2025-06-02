function d = mahal_squared_distances(X, C_inv)  % actually it's squared mahal distance
% Computes SQUARED Mahalanobis distances for each row in X given C_inv
    d = sum((X * C_inv) .* X, 2);  % (X C^{-1}) .* X row-wise sum
end
