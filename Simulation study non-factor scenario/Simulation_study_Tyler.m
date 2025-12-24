%% Run_Tyler_Sim.m
clearvars;
addpath('./GLASD/');           
addpath('./Simulation data/');
addpath('./supp/');
addpath('./other methods/');

num_dataset_rep = 10;         
rng(1)
% ------------------- User Inputs -------------------
p        = 100;
n        = 500;
Ctype    = 'random';    % 'random', 'BlockToeplitz', etc.
dist_type = 't';               % 'gaussian_asymmetric_row', 'gaussian_asymmetric_col', 
                               % 'gaussian_asymmetric_random', 't', etc.
tol_tyler   = 1e-6;
maxit_tyler = 1000;
% ---------------------------------------------------

% Load true correlation matrix
folder = 'Simulation data';
C_true_filename = sprintf('C_p_%d_n_%d_C_%s.csv', p, n, Ctype);
C_true_filepath = fullfile(folder, C_true_filename);

if isfile(C_true_filepath)
    C_true = readmatrix(C_true_filepath);
    fprintf('Loaded true correlation matrix from: %s\n', C_true_filepath);
else
    error('True correlation matrix file not found: %s', C_true_filepath);
end

% Ensure output folder exists
output_folder = 'Simulation output';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

for data_rep = 1:num_dataset_rep
     
    % ---------- Load data matrix X ----------
    folder = 'Simulation data';
    filename = sprintf('X_p_%d_n_%d_C_%s_dist_%s_DataRep_%d.csv', ...
        p, n, Ctype, dist_type, data_rep);
    filepath = fullfile(folder, filename);

    if isfile(filepath)
        X_raw = readmatrix(filepath);
        fprintf('Loaded data from: %s\n', filepath);
    else
        error('File not found: %s', filepath);
    end

    % Use the same standardization routine as in your GLASD pipeline
    X = standardize_data(X_raw);

    % ---------- Compute Tyler's M-estimator ----------
    fprintf('Dataset %d: computing Tyler M-estimator...\n', data_rep);
    C_tyler = Tyler_corr(X, tol_tyler, maxit_tyler);

    % ---------- Compute RMSE and MAD against C_true ----------
    diff_vec = C_tyler(:) - C_true(:);
    rmse_tyler = sqrt(mean(diff_vec.^2));
    mad_tyler  = mean(abs(diff_vec));

    fprintf('Tyler: RMSE = %.4f | MAD = %.4f\n', rmse_tyler, mad_tyler);

    % ---------- Save to CSV in a structure similar to GLASD output ----------
    method_names = {'Tyler'}';   % keep same column style
    T = table(method_names, rmse_tyler, mad_tyler, ...
        'VariableNames', {'Method', 'Best_RMSE', 'Best_MAD'});

    out_filename = fullfile(output_folder, ...
        sprintf('Tyler_rmse_mad_p_%d_n_%d_C_%s_Dist_%s_DataRep_%d.csv', ...
        p, n, Ctype, dist_type, data_rep));

    writetable(T, out_filename);
    fprintf('Saved Tyler results to: %s\n\n', out_filename);

end

fprintf('All Tyler simulations completed for p=%d, n=%d, Ctype=%s, dist=%s.\n', ...
    p, n, Ctype, dist_type);
