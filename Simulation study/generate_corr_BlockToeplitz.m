function C = generate_corr_BlockToeplitz(p, rho_vec, part_prob_vec)
% GENERATE_CORR_BLOCKTOEPLITZ Generates a block Toeplitz correlation matrix.
%
% Inputs:
%   p              - Total dimension of the correlation matrix
%   rho_vec        - Vector of correlation values for each block
%   part_prob_vec  - Vector of probabilities (summing to 1) determining block sizes
%
% Output:
%   C              - p x p block Toeplitz correlation matrix

    if length(rho_vec) ~= length(part_prob_vec)
        error('rho_vec and part_prob_vec must have the same length.');
    end

    if abs(sum(part_prob_vec) - 1) > 1e-6
        error('part_prob_vec must sum to 1.');
    end

    k = length(rho_vec);
    block_sizes = round(part_prob_vec * p);
    diff = p - sum(block_sizes);
    [~, idx] = max(block_sizes);
    block_sizes(idx) = block_sizes(idx) + diff;

    C = eye(p);
    start_idx = 1;

    for b = 1:k
        block_size = block_sizes(b);
        rho = rho_vec(b);
        end_idx = start_idx + block_size - 1;

        for i = start_idx:end_idx
            for j = start_idx:end_idx
                C(i,j) = rho^abs(i-j);
            end
        end

        start_idx = end_idx + 1;
    end
end
