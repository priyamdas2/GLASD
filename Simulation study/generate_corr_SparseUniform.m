function C = generate_corr_SparseUniform(p, sparsity)
% GENERATE_CORR_SPARSEUNIFORM generates a sparse random correlation matrix.
%   p         : number of variables (dimension of the matrix)
%   sparsity  : proportion of off-diagonal entries set to zero (e.g., 0.9)
%
% The nonzero off-diagonal entries are sampled from Uniform[0.1, 0.3].
% All diagonal entries are fixed at 1.

    % Start with identity matrix
    C = eye(p);
    
    % Total number of off-diagonal elements in upper triangle
    num_offdiag = p * (p - 1) / 2;
    num_nonzero = round((1 - sparsity) * num_offdiag);

    % Get linear indices of upper triangle (excluding diagonal)
    upper_inds = find(triu(ones(p), 1));
    selected_inds = randsample(upper_inds, num_nonzero);

    % Generate values from Uniform[0.1, 0.3]
    values = 0.1 + (0.3 - 0.1) * rand(num_nonzero, 1);

    % Fill symmetric off-diagonal values
    for k = 1:num_nonzero
        [i, j] = ind2sub([p, p], selected_inds(k));
        C(i, j) = values(k);
        C(j, i) = values(k);
    end

    % Optional: Project to nearest positive definite matrix
    % (Preserve diagonals at 1)
    [V, D] = eig((C + C') / 2);
    D = max(D, 1e-6);  % enforce positive eigenvalues
    C = V * D * V';
    C = (C + C') / 2;
    C = C ./ sqrt(diag(C) * diag(C)');  % re-standardize to unit diagonals
end
