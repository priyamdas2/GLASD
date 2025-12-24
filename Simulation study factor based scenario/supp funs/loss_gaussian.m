function val = loss_gaussian(X, C)
    C_inv = inv(C);
    d = mahal_squared_distances(X, C_inv);
    val = 0.5 * size(X,1) * log(det(C)) + 0.5 * sum(d);
end
