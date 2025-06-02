clearvars;
addpath('./GLASD/');
addpath('./Simulation data/');
num_dataset_rep = 10;
num_optmz_reps = 10;
params.M       = 100;
params.epsilon = 1e-6;
rng(1)
% ------------------- User Inputs -------------------
p = 20;
n = 100;
Ctype = 'BlockToeplitz';   %  'SparseUniform', 'random', 'BlockToeplitz' (not needed 'Toeplitz', )
dist_type = 't';       % 'gaussian_asymmetric_row', 'gaussian_asymmetric_col', 'gaussian_asymmetric_random', 't'.
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

parfor data_rep = 1:num_dataset_rep
    
    folder = 'Simulation data';
    filename = sprintf('X_p_%d_n_%d_C_%s_dist_%s_DataRep_%d.csv', p, n, Ctype, dist_type, data_rep);
    filepath = fullfile(folder, filename);
    if isfile(filepath)
        X_raw = readmatrix(filepath);
        fprintf('Loaded data from: %s\n', filepath);
    else
        error('File not found: %s', filepath);
    end
    X = standardize_data(X_raw);
    
    
    % Loss 1: Gaussian likelihood
    loss_gaussian_here = @(C) loss_gaussian(X,C);
    
    % Loss 2: Huber
    loss_huber_here = @(C) loss_huber(X,C);
    
    % Loss 3: Truncated
    loss_truncated_here = @(C) loss_truncated(X,C);
    
    % Loss 4: Truncated soft
    loss_truncated_soft_here = @(C) loss_truncated_soft(X,C);
    
    % Loss 5: Tukey Biweight
    loss_tukey_here = @(C) loss_tukey(X,C);
    
    % ---------------------------------------------------
    % Initial point: identity matrix (can be sample corr too)
    C0 = nearestSPD(corr(X));
    d_initial = mahal_squared_distances(X,C0);
    
    % Run GLASD for each loss
    
    results = nan(num_optmz_reps, 10);
    
    % Initialize storage
    obj_vals = nan(num_optmz_reps, 5);  % columns: gauss, huber, trunc, trunc_soft
    C_hat_all = struct('gauss', cell(num_optmz_reps, 1), ...
        'huber', cell(num_optmz_reps, 1), ...
        'trunc', cell(num_optmz_reps, 1), ...
        'trunc_soft', cell(num_optmz_reps, 1),...
        'tukey', cell(num_optmz_reps, 1));
    
    for rep = 1:num_optmz_reps
        rng(rep)
        fprintf('Performing Dataset: %d, iteration: %d\n', data_rep, rep);
        
        C0 = randCorrMatrix(p, rep);
        warning('off', 'all');
        
        C_hat_gauss = GLASD_PD(loss_gaussian_here, C0, params);
        C_hat_huber = GLASD_PD(loss_huber_here, C0, params);
        C_hat_trunc = GLASD_PD(loss_truncated_here, C0, params);
        C_hat_trunc_soft = GLASD_PD(loss_truncated_soft_here, C0,params);
        C_hat_tukey = GLASD_PD(loss_tukey_here, C0,params);
        
        % Store matrices
        C_hat_all(rep).gauss = C_hat_gauss;
        C_hat_all(rep).huber = C_hat_huber;
        C_hat_all(rep).trunc = C_hat_trunc;
        C_hat_all(rep).trunc_soft = C_hat_trunc_soft;
        C_hat_all(rep).tukey = C_hat_tukey;
        
        % Objective values (all evaluated under Gaussian loss)
        obj_vals(rep, 1) = loss_gaussian_here(C_hat_gauss);
        obj_vals(rep, 2) = loss_huber_here(C_hat_huber);
        obj_vals(rep, 3) = loss_truncated_here(C_hat_trunc);
        obj_vals(rep, 4) = loss_truncated_soft_here(C_hat_trunc_soft);
        obj_vals(rep, 5) = loss_tukey_here(C_hat_tukey);
        
    end
    
    % Post-loop: Find the best C_hat (based on lowest obj value for each method)
    [min_vals, min_inds] = min(obj_vals);  % 1 row, 5 columns
    
    best_rmse = zeros(1,5);
    best_mad = zeros(1,5);
    methods = {'gauss', 'huber', 'trunc', 'trunc_soft', 'tukey'};
    
    for i = 1:5
        C_best = C_hat_all(min_inds(i)).(methods{i});
        best_rmse(i) = sqrt(mean((C_best(:) - C_true(:)).^2));
        best_mad(i) = mean(abs(C_best(:) - C_true(:)));
    end
    
    fprintf('\nType of C: %s', Ctype);
    % Display
    fprintf('\nBest Objective Values and Corresponding Errors:\n');
    for i = 1:5
        fprintf('%s: Obj = %.4f | RMSE = %.4f | MAD = %.4f (rep = %d)\n', ...
            methods{i}, min_vals(i), best_rmse(i), best_mad(i), min_inds(i));
    end
    
    method_names = {'gauss', 'huber', 'trunc', 'trunc_soft', 'tukey'}';
    T = table(method_names, best_rmse(:), best_mad(:), ...
        'VariableNames', {'Method', 'Best_RMSE', 'Best_MAD'});
    
    output_folder = 'Simulation output';
    filename = fullfile(output_folder,sprintf('Best_rmse_mad_p_%d_n_%d_NumOptmzReps_%d_C_%s_Dist_%s_DataRep_%d.csv', ...
        p, n, num_optmz_reps, Ctype, dist_type, data_rep));
    writetable(T, filename);
    
end