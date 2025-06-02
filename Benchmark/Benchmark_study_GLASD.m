clearvars; clc; close all;
addpath('./GLASD/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
skip_activeset = 0;
skip_interiorpoint = 0;
skip_sqp = 0;

fun_choices = {'ackley', 'griewank', 'rosenbrock', 'rastrigin'};  % index 1, 2, 3, 4
input_vals = [4];  % example
M = 3;         % M = 5, 10, 20, 50, 100
Num_exp = 2;  % Num_exp = 10 for M = 5, 10, 20, 50, 100.
maxtime_each = 3600; % dont change

GLASD_comp_times = nan(Num_exp,1);
GLASD_funvals = nan(Num_exp,1);

for i = 1:length(input_vals)
    which_fun = fun_choices{input_vals(i)};
    
    if strcmp(which_fun, 'ackley')
        objFun = @(C) modified_ackley(C);
    elseif strcmp(which_fun, 'griewank')
        objFun = @(C) modified_griewank(C);
    elseif strcmp(which_fun, 'rosenbrock')
        objFun = @(C) modified_rosenbrock(C);
    elseif strcmp(which_fun, 'rastrigin')
        objFun = @(C) modified_rastrigin(C);
    else
        error('Unknown function specified in which_fun.');
    end
    
    
    for ii = 1:Num_exp
        fprintf('Performing experiment number: %d using function: %s\n', ii, which_fun);
        rand_seed = ii;
        C0 = randCorrMatrix(M, rand_seed);
        Theta0 = Corr2Theta(C0);
        
        %%%% GLASD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic;
        [C_GLASD, f_GLASD, history] = GLASD_PD(objFun, C0);
        comp_time_GLASD = toc;
        
        
        %%% Summary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        GLASD_comp_times(ii) = comp_time_GLASD;
        
        GLASD_funvals(ii) = f_GLASD;
        
    end
    
    
    % timestamp = datestr(now, 'yyyymmdd_HHMM');
    
    output_folder = 'ALL_Outputs';
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end
    
    filename = fullfile(output_folder, ['Comp_times_GLASD_' which_fun '_M_' num2str(M) '_reps_' num2str(Num_exp) '.csv']);
    writematrix(GLASD_comp_times, filename);
    
    filename = fullfile(output_folder, ['Funvals_GLASD_' which_fun '_M_' num2str(M) '_reps_' num2str(Num_exp) '.csv']);
    writematrix(GLASD_funvals, filename);
    
    
end

