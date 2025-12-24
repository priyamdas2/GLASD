clearvars;
addpath('./GLASD/');
addpath('./Real data/');
addpath('./supp/');
addpath('./other methods/');
rng(1)
X_raw = readmatrix('X_blocked_order.csv');


%% Naive estimate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_naive = corr(X_raw);
writematrix(C_naive, 'C_naive.csv');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Tyler estimate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol_tyler   = 1e-6;
maxit_tyler = 1000;

C_tyler = Tyler_corr(X_raw, tol_tyler, maxit_tyler);
writematrix(C_tyler, 'C_tyler.csv');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% POET estimate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts = struct();
opts.center    = true;
opts.tau_const = 2;        % same universal threshold as FWZ
opts.threshold = 'hard';  % hard thresholding (standard POET)

r = select_r_BaiNg(X_raw, 100);
fit1 = poet2013_cov_estimator(X_raw, r, opts);
Cov_POET   = (fit1.Sigma_hat   + fit1.Sigma_hat')   / 2;

diag_mat = sqrt(diag(Cov_POET));
C_POET = Cov_POET ./ (diag_mat * diag_mat');
writematrix(C_POET, 'C_POET.csv');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


X = standardize_data(X_raw);
p = size(X,2);
num_reps = 1;

loss_huber_here = @(C) loss_huber(X, C);

% Store loss values and corresponding solutions
all_losses = zeros(num_reps, 1);
all_solutions = cell(num_reps, 1);
Time_taken = zeros(num_reps, 1);
for rep = 1:num_reps
    rng(rep)
    tic;
    fprintf('Running repetition %d of %d...\n', rep, num_reps);
    C0 = nearest_corr_fast(C_POET);  % initial point
    warning('off', 'all');
    
    params.M       = 100;
    C_hat_huber = GLASD_PD(loss_huber_here, C0, params);
    loss_val = loss_huber_here(C_hat_huber);  % evaluate loss
    
    all_losses(rep) = loss_val;
    all_solutions{rep} = C_hat_huber;
    Time_taken(rep) = toc;
end

% Select best solution
[~, best_idx] = min(all_losses);
C_optimal = all_solutions{best_idx};

writematrix(C_optimal, 'C_huber_optimal.csv');
writematrix(Time_taken, 'Time_taken_GLASD.csv');