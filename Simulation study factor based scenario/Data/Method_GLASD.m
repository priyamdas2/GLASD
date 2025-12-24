%% ============================================
% GLASD Model Fitting & Errors
% Fits GLASD (5 losses) from Y only
% Computes Σ / Σ^{-1} errors vs true Σ
% SINGLE-START (no multistart)
% ============================================

clear; clc;

addpath('./GLASD/');
addpath('./supp funs/');

%% --------------------------------------------
% Simulation identifiers (must match data script)
% --------------------------------------------

p   = 50;
r   = 3;
n   = p / 2;
nu  = 3;
num_exp = 10;
Sigma_u_regime = 'random';   % 'sparse' / 'diagonal' / 'random'

for rep = 1:num_exp
    %% --------------------------------------------
    % Load B and Sigma_u (to reconstruct true Sigma)
    % --------------------------------------------
    
    fname_B  = sprintf('Data/B_matrix_p_%d_nu_%d_rep_%d.csv', p, nu, rep);
    B        = readmatrix(fname_B);
    
    fname_Su = sprintf('Data/SigmaU_%s_p_%d_nu_%d_rep_%d.csv', Sigma_u_regime, p, nu, rep);
    Sigma_u_true = readmatrix(fname_Su);
    
    Sigma_true     = B * B' + Sigma_u_true;
    Sigma_inv_true = inv(Sigma_true);
    
    methods = {'gauss', 'huber', 'trunc', 'trunc_soft', 'tukey'};
    
    %% ============================================
    % Setting 1: Elliptical multivariate t
    % ============================================
    
    fname_Y1 = sprintf('Data/Y_%s_setting1_elliptical_p_%d_nu_%d_rep_%d.csv', ...
        Sigma_u_regime, p, nu, rep);
    Y1 = readmatrix(fname_Y1);   % n x p
    
    % Standardize data as in GLASD sims
    X1 = standardize_data(Y1);
    
    % Estimated marginal scales (NOT oracle)
    D_hat_1 = diag(std(Y1));
    
    % Define loss functions for X1
    loss_gaussian_1       = @(C) loss_gaussian(X1, C);
    loss_huber_1          = @(C) loss_huber(X1, C);
    loss_truncated_1      = @(C) loss_truncated(X1, C);
    loss_truncated_soft_1 = @(C) loss_truncated_soft(X1, C);
    loss_tukey_1          = @(C) loss_tukey(X1, C);
    
    loss_list_1 = {loss_gaussian_1, loss_huber_1, loss_truncated_1, ...
        loss_truncated_soft_1, loss_tukey_1};
    
    % Error matrix: 5 methods x 4 metrics
    % [||Σ̂-Σ||_F, ||Σ̂-Σ||_2, ||Σ̂^{-1}-Σ^{-1}||_F, ||Σ̂^{-1}-Σ^{-1}||_2]
    err_glasd_1 = nan(5, 4);
    
    for m = 1:5
        fprintf('Performing experiment number: %d, Setting: 1, Loss func: %d.\n', rep, m);
        if(m > 1)
           fprintf('Last loss func in this scenario took %.2f secs.\n', time_taken);
        end
        tic;
        loss_fun = loss_list_1{m};
        
        % Single deterministic start
        rng(1);
        C0 = randCorrMatrix(p, 1);
        warning('off', 'all');
        
        C_hat = GLASD_PD(loss_fun, C0);
        
        % Symmetrize correlation estimate
        C_best = (C_hat + C_hat') / 2;
        
        % Map corr -> cov using ESTIMATED marginal variances
        Sigma_hat = D_hat_1 * C_best * D_hat_1;
        Sigma_hat = (Sigma_hat + Sigma_hat') / 2;
        
        % Invert WITHOUT artificial regularization
        Sigma_inv_hat = inv(Sigma_hat);
        
        % Errors vs true Σ and Σ^{-1}
        err_Sigma_frob = norm(Sigma_hat     - Sigma_true,     'fro');
        err_Sigma_op   = norm(Sigma_hat     - Sigma_true);
        err_Prec_frob  = norm(Sigma_inv_hat - Sigma_inv_true, 'fro');
        err_Prec_op    = norm(Sigma_inv_hat - Sigma_inv_true);
        
        err_glasd_1(m, :) = [err_Sigma_frob, err_Sigma_op, err_Prec_frob, err_Prec_op];
        time_taken = toc;
    end
    
    fname_err1 = sprintf('Output/GLASD_errors_%s_setting1_elliptical_p_%d_nu_%d_rep_%d.csv', ...
        Sigma_u_regime, p, nu, rep);
    writematrix(err_glasd_1, fname_err1);
    
    %% ============================================
    % Setting 2: Componentwise iid t
    % ============================================
    
    fname_Y2 = sprintf('Data/Y_%s_setting2_componentwise_p_%d_nu_%d_rep_%d.csv', ...
        Sigma_u_regime, p, nu, rep);
    Y2 = readmatrix(fname_Y2);   % n x p
    
    X2 = standardize_data(Y2);
    
    D_hat_2 = diag(std(Y2));
    
    loss_gaussian_2       = @(C) loss_gaussian(X2, C);
    loss_huber_2          = @(C) loss_huber(X2, C);
    loss_truncated_2      = @(C) loss_truncated(X2, C);
    loss_truncated_soft_2 = @(C) loss_truncated_soft(X2, C);
    loss_tukey_2          = @(C) loss_tukey(X2, C);
    
    loss_list_2 = {loss_gaussian_2, loss_huber_2, loss_truncated_2, ...
        loss_truncated_soft_2, loss_tukey_2};
    
    err_glasd_2 = nan(5, 4);
    
    for m = 1:5
        fprintf('Performing experiment number: %d, Setting: 2, Loss func: %d.\n', rep, m);
        if(m > 1)
           fprintf('Last loss func in this scenario took %.2f secs.\n', time_taken);
        end
        
        loss_fun = loss_list_2{m};
        
        rng(1);
        C0 = randCorrMatrix(p, 1);
        warning('off', 'all');
        
        C_hat = GLASD_PD(loss_fun, C0);
        
        C_best = (C_hat + C_hat') / 2;
        
        Sigma_hat = D_hat_2 * C_best * D_hat_2;
        Sigma_hat = (Sigma_hat + Sigma_hat') / 2;
        
        Sigma_inv_hat = inv(Sigma_hat);
        
        err_Sigma_frob = norm(Sigma_hat     - Sigma_true,     'fro');
        err_Sigma_op   = norm(Sigma_hat     - Sigma_true);
        err_Prec_frob  = norm(Sigma_inv_hat - Sigma_inv_true, 'fro');
        err_Prec_op    = norm(Sigma_inv_hat - Sigma_inv_true);
        
        err_glasd_2(m, :) = [err_Sigma_frob, err_Sigma_op, err_Prec_frob, err_Prec_op];
    end
    
    fname_err2 = sprintf('Output/GLASD_errors_%s_setting2_componentwise_p_%d_nu_%d_rep_%d.csv', ...
        Sigma_u_regime, p, nu, rep);
    writematrix(err_glasd_2, fname_err2);
    
    %% ============================================
    % Console confirmation
    % ============================================
    
    disp('GLASD error files written:');
    disp(['  ', fname_err1]);
    disp(['  ', fname_err2]);
    
end
