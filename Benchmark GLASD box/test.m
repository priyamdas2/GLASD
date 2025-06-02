clear all
rng(3)
n = 100;

% % Custom GLASD parameters for Ackley
% params.s_init  = 0.5;
% params.p_init  = ones(2*n,1);  % Uniform initial prob
% params.s_inc   = 2;
% params.s_dec   = 2;
% params.p_inc   = 2;
% params.p_dec   = 2;
% params.m       = 5;              % More exploration
% params.c       = 0.0001*log(n);         % Cooling constant
% params.T       = 3000*log(n);    % Max iterations
% params.M       = 4*n;            % Stagnation window
% params.epsilon = 1e-20;           % Convergence tolerance


%%% Ackley %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lb = -32.768*ones(n,1);
ub = 32.768*ones(n,1);
x0 = lb + (ub - lb) .* rand(n, 1); 

% Run GLASD
[x_opt, f_opt, hist] = GLASD(@ackley, lb, ub, x0, params);

% Display results
fprintf('Best function value: %d\n', f_opt);
disp('Best solution found:');
%disp(x_opt');

epsilon = 1e-20;  % to avoid log(0) errors
log_fvals = log10(hist.fvals + epsilon);

figure;
plot(log_fvals, 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('log_{10}(Best function value)');
title('GLASD on Ackley Function');
grid on;


%%% Griewank %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lb = -600*ones(n,1);
ub = 600*ones(n,1);
x0 = lb + (ub - lb) .* rand(n, 1); 

% Run GLASD
[x_opt, f_opt, hist] = GLASD(@griewank, lb, ub, x0, params);

% Display results
fprintf('Best function value: %d\n', f_opt);
disp('Best solution found:');
%disp(x_opt');

epsilon = 1e-20;  % to avoid log(0) errors
log_fvals = log10(hist.fvals + epsilon);

figure;
plot(log_fvals, 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('log_{10}(Best function value)');
title('GLASD on Griewank Function');
grid on;

%%% Rastrigin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lb = -5.12*ones(n,1);
ub = 5.12*ones(n,1);
x0 = lb + (ub - lb) .* rand(n, 1); 

% Run GLASD
[x_opt, f_opt, hist] = GLASD(@rastrigin, lb, ub, x0, params);

% Display results
fprintf('Best function value: %d\n', f_opt);
disp('Best solution found:');
%disp(x_opt');

epsilon = 1e-20;  % to avoid log(0) errors
log_fvals = log10(hist.fvals + epsilon);

figure;
plot(log_fvals, 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('log_{10}(Best function value)');
title('GLASD on Rastrigin Function');
grid on;