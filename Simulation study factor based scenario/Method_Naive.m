%% ============================================
% Naive Sample Covariance Model Fitting & Errors
% Baseline: Σ̂ = sample covariance of Y
% Computes Σ / Σ^{-1} errors vs true Σ
%
% OUTPUT (2 CSVs):
%   Data/Naive_errors_<regime>_setting1_elliptical_p_*_nu_*_rep_*.csv
%   Data/Naive_errors_<regime>_setting2_componentwise_p_*_nu_*_rep_*.csv
% ============================================

clear; clc;

addpath('./supp funs/');

%% --------------------------------------------
% Simulation identifiers (must match data script)
% --------------------------------------------

p   = 50;
r   = 3;
n   = p / 2;
nu  = 3;
num_exp = 10;
Sigma_u_regime = 'random';   % 'diagonal' / 'random'

for rep = 1:num_exp    
    %% --------------------------------------------
    % Load B and Sigma_u (to reconstruct true Sigma)
    % --------------------------------------------
    
    fname_B  = sprintf('Data/B_matrix_p_%d_nu_%d_rep_%d.csv', p, nu, rep);
    B        = readmatrix(fname_B);
    
    fname_Su = sprintf('Data/SigmaU_%s_p_%d_nu_%d_rep_%d.csv', ...
        Sigma_u_regime, p, nu, rep);
    Sigma_u_true = readmatrix(fname_Su);
    
    Sigma_true     = B * B' + Sigma_u_true;
    Sigma_inv_true = inv(Sigma_true);
    
    %% ============================================
    % Setting 1: Elliptical multivariate t
    % ============================================
    
    fname_Y1 = sprintf('Data/Y_%s_setting1_elliptical_p_%d_nu_%d_rep_%d.csv', ...
        Sigma_u_regime, p, nu, rep);
    Y1 = readmatrix(fname_Y1);   % n x p
    
    % --- Naive estimator ---
    Sigma_hat_1 = cov(Y1, 1);
    Sigma_hat_1 = (Sigma_hat_1 + Sigma_hat_1') / 2;
    
    Sigma_inv_hat_1 = inv(Sigma_hat_1);
    
    % --- Error metrics ---
    err_Sigma_frob_1 = norm(Sigma_hat_1     - Sigma_true,     'fro');
    err_Sigma_op_1   = norm(Sigma_hat_1     - Sigma_true);
    err_Prec_frob_1  = norm(Sigma_inv_hat_1 - Sigma_inv_true, 'fro');
    err_Prec_op_1    = norm(Sigma_inv_hat_1 - Sigma_inv_true);
    
    err_naive_1 = [err_Sigma_frob_1, err_Sigma_op_1, ...
        err_Prec_frob_1,  err_Prec_op_1];
    
    fname_err1 = sprintf(['Output/Naive_errors_%s_setting1_elliptical_' ...
        'p_%d_nu_%d_rep_%d.csv'], ...
        Sigma_u_regime, p, nu, rep);
    writematrix(err_naive_1, fname_err1);
    
    
    %% ============================================
    % Setting 2: Componentwise iid t
    % ============================================
    
    fname_Y2 = sprintf('Data/Y_%s_setting2_componentwise_p_%d_nu_%d_rep_%d.csv', ...
        Sigma_u_regime, p, nu, rep);
    Y2 = readmatrix(fname_Y2);   % n x p
    
    Sigma_hat_2 = cov(Y2, 1);
    Sigma_hat_2 = (Sigma_hat_2 + Sigma_hat_2') / 2;
    
    Sigma_inv_hat_2 = inv(Sigma_hat_2);
    
    err_Sigma_frob_2 = norm(Sigma_hat_2     - Sigma_true,     'fro');
    err_Sigma_op_2   = norm(Sigma_hat_2     - Sigma_true);
    err_Prec_frob_2  = norm(Sigma_inv_hat_2 - Sigma_inv_true, 'fro');
    err_Prec_op_2    = norm(Sigma_inv_hat_2 - Sigma_inv_true);
    
    err_naive_2 = [err_Sigma_frob_2, err_Sigma_op_2, ...
        err_Prec_frob_2,  err_Prec_op_2];
    
    fname_err2 = sprintf(['Output/Naive_errors_%s_setting2_componentwise_' ...
        'p_%d_nu_%d_rep_%d.csv'], ...
        Sigma_u_regime, p, nu, rep);
    writematrix(err_naive_2, fname_err2);
    
    
    %% ============================================
    % Console confirmation
    % ============================================
    
    disp('Naive covariance error files written:');
    disp(['  ', fname_err1]);
    disp(['  ', fname_err2]);
    
end
