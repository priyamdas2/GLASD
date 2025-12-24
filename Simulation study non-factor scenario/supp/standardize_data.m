function X_std = standardize_data(X)
% STANDARDIZE_DATA standardizes the columns of X to have zero mean and unit variance.
%
%   Input:
%     X      : n x p data matrix (n observations, p variables)
%
%   Output:
%     X_std  : standardized data matrix

    mu = mean(X);           % 1 x p vector of column means
    sigma = std(X);         % 1 x p vector of column std deviations
    X_std = (X - mu) ./ sigma;  % elementwise standardization
end
