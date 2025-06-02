clearvars;
addpath('./GLASD/');
addpath('./Simulation data/');

X_raw = readmatrix('BRCA_subset.csv');
X = standardize_data(X_raw);
p = size(X,2);
num_reps = 50;

loss_huber_here = @(C) loss_huber(X, C);

% Store loss values and corresponding solutions
all_losses = zeros(num_reps, 1);
all_solutions = cell(num_reps, 1);

for rep = 1:num_reps
    fprintf('Running repetition %d of %d...\n', rep, num_reps);
    C0 = randCorrMatrix(p, rep);  % initial point
    warning('off', 'all');
    
    params.M       = 100;
    params.epsilon = 1e-6;
    C_hat_huber = GLASD_PD(loss_huber_here, C0, params);
    loss_val = loss_huber_here(C_hat_huber);  % evaluate loss
    
    all_losses(rep) = loss_val;
    all_solutions{rep} = C_hat_huber;
end

% Select best solution
[~, best_idx] = min(all_losses);
C_optimal = all_solutions{best_idx};

writematrix(C_optimal, 'C_huber_optimal.csv');