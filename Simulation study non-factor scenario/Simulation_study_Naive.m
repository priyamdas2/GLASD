%% Run_Naive_Sim.m
clearvars;
addpath('./GLASD/');           
addpath('./Simulation data/');
addpath('./supp/');
addpath('./other methods/');   % harmless if you keep everything consistent

num_dataset_rep = 10;         
rng(1)

% ------------------- User Inputs -------------------
p         = 100;
n         = 500;
Ctype     = 'random';    % 'random', 'BlockToeplitz', etc.
dist_type = 'gaussian_asymmetric_col';  % 'gaussian_asymmetric_row', 
                                           % 'gaussian_asymmetric_col', 
                                           % 'gaussian_asymmetric_random', 't', etc.
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

    % ---------- Compute Naive sample correlation ----------
    fprintf('Dataset %d: computing Naive sample correlation...\n', data_rep);

    % Sample correlation (p x p)
    C_naive = corr(X);

    % (Optional, for numerical safety you could project to SPD)
    % if exist('nearestSPD', 'file') == 2
    %     C_naive = nearestSPD(C_naive);
    % end

    % ---------- Compute RMSE and MAD against C_true ----------
    diff_vec   = C_naive(:) - C_true(:);
    rmse_naive = sqrt(mean(diff_vec.^2));
    mad_naive  = mean(abs(diff_vec));

    fprintf('Naive: RMSE = %.4f | MAD = %.4f\n', rmse_naive, mad_naive);

    % ---------- Save to CSV in a structure similar to GLASD/Tyler ----------
    method_names = {'Naive'}';   % keep same column style
    T = table(method_names, rmse_naive, mad_naive, ...
        'VariableNames', {'Method', 'Best_RMSE', 'Best_MAD'});

    out_filename = fullfile(output_folder, ...
        sprintf('Naive_rmse_mad_p_%d_n_%d_C_%s_Dist_%s_DataRep_%d.csv', ...
        p, n, Ctype, dist_type, data_rep));

    writetable(T, out_filename);
    fprintf('Saved Naive results to: %s\n\n', out_filename);

end

fprintf('All Naive simulations completed for p=%d, n=%d, Ctype=%s, dist=%s.\n', ...
    p, n, Ctype, dist_type);
