rng(1)

% -------------------- USER SETTINGS --------------------
p = 50;                      % number of variables
n = 500;                      % number of samples
num_data_rep = 10;
structure_options = [2 3 4];  % choose from: 1 = Toeplitz, 2 = BlockToeplitz, 3 = SparseUniform, 4 = Random
dist_type = 'gaussian_asymmetric_random';      % choose: % 'gaussian_asymmetric_row', 'gaussian_asymmetric_col', 'gaussian_asymmetric_random', 't'.

% -------------------------------------------------------

% Create output folder if it doesn't exist
output_folder = 'Simulation data';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

for i = 1:length(structure_options)
    opt = structure_options(i);
    switch opt
        case 1
            C = generate_corr_Toeplitz(p, 0.6);
            name = 'Toeplitz';
        case 2
            rho_vec = [0.6, 0.3, 0.4];
            part_prob_vec = [0.25, 0.5, 0.25];
            C = generate_corr_BlockToeplitz(p, rho_vec, part_prob_vec);
            name = 'BlockToeplitz';
        case 3
            C = generate_corr_SparseUniform(p, 0.9);
            name = 'SparseUniform';
        case 4
            C = randCorrMatrix(p, 123);
            name = 'random';
        otherwise
            error('Unknown structure option: %d', opt);
    end
    filename_C = fullfile(output_folder, sprintf('C_p_%d_n_%d_C_%s.csv', p, n, name));
    writematrix(C, filename_C);
    
    for data_rep = 1:num_data_rep
        rng(data_rep)
        % Simulate data
        X = simulate_data_from_corr(C, n, dist_type);

        filename_X = fullfile(output_folder, sprintf('X_p_%d_n_%d_C_%s_dist_%s_DataRep_%d.csv', p, n, name, dist_type, data_rep));
        writematrix(X, filename_X);
    end
    
    
end
