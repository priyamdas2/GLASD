function C = generate_corr_Toeplitz(p, rho)
% GEN_CORR_TOEPLITZ Generates a p x p Toeplitz correlation matrix
% with entries C_ij = rho^|i - j|.
%
% Inputs:
%   rho : scalar in (0, 1), correlation decay parameter
%   p   : integer, dimension of the correlation matrix
%
% Output:
%   C   : p x p Toeplitz correlation matrix

    if rho <= 0 || rho >= 1
        error('rho must be in the open interval (0, 1).');
    end

    if p <= 1 || floor(p) ~= p
        error('p must be an integer greater than 1.');
    end

    % Generate first row of the Toeplitz matrix
    first_row = rho .^ (0:(p - 1));

    % Construct the Toeplitz matrix
    C = toeplitz(first_row);

end