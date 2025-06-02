clear; clc; close all;
rng(1)
n = 100;
lb = -5 * ones(n,1);
ub = 10 * ones(n,1);
x0 = lb + (ub - lb) .* rand(n, 1);

% Setup global eval counter
global eval_count eval_history
eval_count = 0;
eval_history = [];

% Define inline rosenbrock with counter
rosenbrock_eval = @(x) eval_count_store(rosenbrock(x));
rosenbrock_eval_ga = @(x) arrayfun(@(i) eval_count_store(rosenbrock(x(i,:)')), 1:size(x,1))';
% GLASD
eval_count = 0;
params.T = round(3000 * log(n));

[xg, fg, histg] = GLASD(rosenbrock_eval, lb, ub, x0, params);
evals_glasd = 0:length(histg.fvals)-1;

% Simulated Annealing
eval_count = 0;
eval_history = [];
options_sa = saoptimset('Display','off','MaxIter',params.T, ...
    'PlotFcns',{}, ...
    'OutputFcn', @(x, optimValues, state) deal(false, optimValues, false));
[xs, fs] = simulannealbnd(rosenbrock_eval, x0, lb, ub, options_sa);
evals_sa = eval_history;

% Genetic Algorithm
eval_count = 0;
eval_history = [];

options_ga = optimoptions('ga', ...
    'Display','off', ...
    'MaxGenerations',500, ...
    'PlotFcn',{}, ...
    'UseParallel', false, ...
    'InitialPopulationMatrix', x0(:)');  

[xga, fga] = ga(rosenbrock_eval_ga, n, [], [], [], [], lb, ub, [], options_ga);
evals_ga = eval_history;

% Particle Swarm Optimization (PSO)

eval_count = 0;
eval_history = [];

options_pso = optimoptions('particleswarm', ...
    'Display', 'off', ...
    'MaxIterations', params.T);

[xpso, fpso] = particleswarm(rosenbrock_eval_ga, n, lb, ub, options_pso);
evals_pso = eval_history;

% Pattern Search
eval_count = 0;
eval_history = [];
options_ps = optimoptions('patternsearch', ...
    'Display', 'off', ...
    'MaxIterations', params.T, ...
    'PlotFcn', {});
[xps, fps] = patternsearch(rosenbrock_eval, x0, [], [], [], [], lb, ub, options_ps);
evals_ps = eval_history;

% Interior-point
eval_count = 0;
eval_history = [];
options_fmc1 = optimoptions('fmincon', ...
    'Display','off', ...
    'Algorithm','interior-point', ...
    'MaxIterations', params.T);
[xfmc1, ffmc1] = fmincon(rosenbrock_eval, x0, [], [], [], [], lb, ub, [], options_fmc1);
evals_fmc1 = eval_history;

% SQP
eval_count = 0;
eval_history = [];
options_fmc2 = optimoptions('fmincon', ...
    'Display','off', ...
    'Algorithm','sqp', ...
    'MaxIterations', params.T);
[xfmc2, ffmc2] = fmincon(rosenbrock_eval, x0, [], [], [], [], lb, ub, [], options_fmc2);
evals_fmc2 = eval_history;

% Active-set
eval_count = 0;
eval_history = [];
options_fmc3 = optimoptions('fmincon', ...
    'Display','off', ...
    'Algorithm','active-set', ...
    'MaxIterations', params.T);
[xfmc3, ffmc3] = fmincon(rosenbrock_eval, x0, [], [], [], [], lb, ub, [], options_fmc3);
evals_fmc3 = eval_history;



% Convert objective values to log scale
eps_val = 1e-20;
plot_glasd = log(histg.fvals + eps_val);
plot_sa    = log(cummin(evals_sa) + eps_val);
plot_ga    = log(cummin(evals_ga) + eps_val);
plot_pso   = log(cummin(evals_pso) + eps_val);
plot_ps    = log(cummin(evals_ps) + eps_val);
plot_fmc1  = log(cummin(evals_fmc1) + eps_val);
plot_fmc2  = log(cummin(evals_fmc2) + eps_val);
plot_fmc3  = log(cummin(evals_fmc3) + eps_val);

%%% Plot without legend %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
hold on;

plot(evals_glasd, plot_glasd, 'LineWidth', 1.2, 'Color', [0, 0, 1]);
plot(1:length(plot_sa),   plot_sa,   'LineWidth', 1.2, 'Color', [0.85, 0.33, 0.1]);
plot(1:length(plot_ga),   plot_ga,   'LineWidth', 1.2, 'Color', [0.47, 0.67, 0.19]);
plot(1:length(plot_pso),  plot_pso,  'LineWidth', 1.2, 'Color', [0.49, 0.18, 0.56]);
plot(1:length(plot_ps),   plot_ps,   'LineWidth', 1.2, 'Color', [0.3, 0.75, 0.93]);
plot(1:length(plot_fmc1), plot_fmc1, 'LineWidth', 1.2, 'Color', [0.93, 0.69, 0.13]);
plot(1:length(plot_fmc2), plot_fmc2, 'LineWidth', 1.2, 'Color', [0.64, 0.08, 0.18]);
plot(1:length(plot_fmc3), plot_fmc3, 'LineWidth', 1.2, 'Color', [0.2, 0.2, 0.2]);


xlabel('Function evaluations', 'FontSize', 12);
ylabel('$\log_e(\mathrm{Objective\ Value})$', 'Interpreter', 'latex', ...
       'FontSize', 14, 'FontWeight', 'bold');
title('Rosenbrock', 'FontSize', 16, 'FontWeight', 'bold');
xlim([0 4000]);
grid on;

% Save plot WITHOUT legend
print(gcf, 'plot_rosenbrock', '-dpng', '-r600');  % 600 DPI PNG


%%% Plot of LEGEND only %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;
% Dummy plots for legend only
h1 = plot(nan, nan, '-', 'LineWidth', 1.8, 'Color', [0, 0, 1]);
h2 = plot(nan, nan, '-', 'LineWidth', 1.8, 'Color', [0.85, 0.33, 0.1]);
h3 = plot(nan, nan, '-', 'LineWidth', 1.8, 'Color', [0.47, 0.67, 0.19]);
h4 = plot(nan, nan, '-', 'LineWidth', 1.8, 'Color', [0.49, 0.18, 0.56]);
h5 = plot(nan, nan, '-', 'LineWidth', 1.8, 'Color', [0.3, 0.75, 0.93]);
h6 = plot(nan, nan, '-', 'LineWidth', 1.8, 'Color', [0.93, 0.69, 0.13]);
h7 = plot(nan, nan, '-', 'LineWidth', 1.8, 'Color', [0.64, 0.08, 0.18]);
h8 = plot(nan, nan, '-', 'LineWidth', 1.8, 'Color', [0.2, 0.2, 0.2]);

legend([h1 h2 h3 h4 h5 h6 h7 h8], ...
    {'GLASD','Simulated Annealing','Genetic Algorithm', ...
     'Particle Swarm Opt','Pattern Search','Interior Point', ...
     'Seq Quad Prog','Active Set'}, ...
     'FontSize', 12, 'Location', 'northwest');

axis off;
set(gcf, 'Position', [400 200 300 300]);  % Adjust size

% Save legend figure
print(gcf, 'plot_legend_only', '-dpng', '-r600');