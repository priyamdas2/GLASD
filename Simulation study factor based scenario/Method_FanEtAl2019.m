%% ============================================
% Fan–Wang–Zhong (2019) Model Fitting & Errors
% Computes error metrics for FWZ estimator only
% (no GLASD fitting here)
% ============================================

clear; clc;

addpath('./FanEtAl2019/');
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
    
    fname_Su = sprintf('Data/SigmaU_%s_p_%d_nu_%d_rep_%d.csv', Sigma_u_regime, p, nu, rep);
    Sigma_u_true = readmatrix(fname_Su);
    
    Sigma_true     = B * B' + Sigma_u_true;
    Sigma_inv_true = inv(Sigma_true);
    
    % Off-diagonal mask for Sigma_u errors
    mask = ~eye(p);
    
    %% --------------------------------------------
    % FWZ tuning parameters
    % --------------------------------------------
    
    opts = struct();
    opts.center          = true;
    opts.tau_const       = 2;
    opts.alpha_mult_diag = 1.0;
    opts.alpha_mult_off  = 0.8;
    
    %% ============================================
    % Setting 1: Elliptical multivariate t
    % ============================================
    
    fname_Y1 = sprintf('Data/Y_%s_setting1_elliptical_p_%d_nu_%d_rep_%d.csv', ...
        Sigma_u_regime, p, nu, rep);
    fname_F1 = sprintf('Data/F_%s_setting1_elliptical_p_%d_nu_%d_rep_%d.csv', ...
        Sigma_u_regime, p, nu, rep);
    
    Y1 = readmatrix(fname_Y1);   % n x p
    F1 = readmatrix(fname_F1);   % n x r
    
    % ---- Fit Fan et al. 2019 (FWZ) on setting 1 ----
    fit1 = fan2019_cov_estimator(Y1, F1, opts);
    
    Sigma_hat_1   = (fit1.Sigma_hat   + fit1.Sigma_hat')   / 2;
    Sigma_u_hat_1 = (fit1.Sigma_u_hat + fit1.Sigma_u_hat') / 2;
    
    % ---- Inversion (no artificial eigen regularization) ----
    Sigma_inv_hat_1 = inv(Sigma_hat_1);
    
    % ---- Error metrics for setting 1 ----
    err_Sigma_frob_1   = norm(Sigma_hat_1     - Sigma_true,     'fro');
    err_Sigma_op_1     = norm(Sigma_hat_1     - Sigma_true);
    err_Prec_frob_1    = norm(Sigma_inv_hat_1 - Sigma_inv_true, 'fro');
    err_Prec_op_1      = norm(Sigma_inv_hat_1 - Sigma_inv_true);
    err_Sigma_u_frob_1 = norm((Sigma_u_hat_1  - Sigma_u_true) .* mask, 'fro');
    err_Sigma_u_op_1   = norm((Sigma_u_hat_1  - Sigma_u_true) .* mask);
    
    err1 = [err_Sigma_frob_1, err_Sigma_op_1, ...
        err_Prec_frob_1,  err_Prec_op_1, ...
        err_Sigma_u_frob_1, err_Sigma_u_op_1];
    
    fname_err1 = sprintf('Output/FWZ_errors_%s_setting1_elliptical_p_%d_nu_%d_rep_%d.csv', ...
        Sigma_u_regime, p, nu, rep);
    writematrix(err1, fname_err1);
    
    %% ============================================
    % Setting 2: Componentwise iid t
    % ============================================
    
    fname_Y2 = sprintf('Data/Y_%s_setting2_componentwise_p_%d_nu_%d_rep_%d.csv', ...
        Sigma_u_regime, p, nu, rep);
    fname_F2 = sprintf('Data/F_%s_setting2_componentwise_p_%d_nu_%d_rep_%d.csv', ...
        Sigma_u_regime, p, nu, rep);
    
    Y2 = readmatrix(fname_Y2);   % n x p
    F2 = readmatrix(fname_F2);   % n x r
    
    % ---- Fit Fan et al. 2019 (FWZ) on setting 2 ----
    fit2 = fan2019_cov_estimator(Y2, F2, opts);
    
    Sigma_hat_2   = (fit2.Sigma_hat   + fit2.Sigma_hat')   / 2;
    Sigma_u_hat_2 = (fit2.Sigma_u_hat + fit2.Sigma_u_hat') / 2;
    
    Sigma_inv_hat_2 = inv(Sigma_hat_2);
    
    % ---- Error metrics for setting 2 ----
    err_Sigma_frob_2   = norm(Sigma_hat_2     - Sigma_true,     'fro');
    err_Sigma_op_2     = norm(Sigma_hat_2     - Sigma_true);
    err_Prec_frob_2    = norm(Sigma_inv_hat_2 - Sigma_inv_true, 'fro');
    err_Prec_op_2      = norm(Sigma_inv_hat_2 - Sigma_inv_true);
    err_Sigma_u_frob_2 = norm((Sigma_u_hat_2  - Sigma_u_true) .* mask, 'fro');
    err_Sigma_u_op_2   = norm((Sigma_u_hat_2  - Sigma_u_true) .* mask);
    
    err2 = [err_Sigma_frob_2, err_Sigma_op_2, ...
        err_Prec_frob_2,  err_Prec_op_2, ...
        err_Sigma_u_frob_2, err_Sigma_u_op_2];
    
    fname_err2 = sprintf('Output/FWZ_errors_%s_setting2_componentwise_p_%d_nu_%d_rep_%d.csv', ...
        Sigma_u_regime, p, nu, rep);
    writematrix(err2, fname_err2);
    
    disp('Robust POET error files written:');
    disp(['  ', fname_err1]);
    disp(['  ', fname_err2]);
end
