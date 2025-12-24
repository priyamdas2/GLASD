% ==== USER INPUTS ====
p = 100;
n = 500;
NumOptmzReps   = 10;
NumDatasetReps = 10;
Ctype = 'random'; % 'random', 'BlockToeplitz'
Dist  = 'gaussian_asymmetric_random';  % 'gaussian_asymmetric_row', 'gaussian_asymmetric_col', 'gaussian_asymmetric_random', 't'.

% ==== FILE READING LOOP FOR GLASD OUTPUTS ====
results    = cell(1, NumDatasetReps);   % to store each matrix
valid_idx  = false(1, NumDatasetReps);  % which dataset reps are actually found

for DataRep = 1:NumDatasetReps
    % Construct filename (GLASD summary for this dataset)
    filename = sprintf('Best_rmse_mad_p_%d_n_%d_NumOptmzReps_%d_C_%s_Dist_%s_DataRep_%d.csv', ...
        p, n, NumOptmzReps, Ctype, Dist, DataRep);
    
    if isfile(filename)
        results{DataRep} = readmatrix(filename);
        valid_idx(DataRep) = true;
    else
        warning('File not found: %s', filename);
        results{DataRep} = [];
    end
end

% ==== FILTER VALID MATRICES ====
valid_results = results(valid_idx);

if isempty(valid_results)
    error('No valid GLASD matrices found to summarize.');
end

% ==== VERIFY ALL GLASD MATRICES SAME SIZE ====
matrix_size = size(valid_results{1});
if ~all(cellfun(@(x) isequal(size(x), matrix_size), valid_results))
    error('Inconsistent matrix sizes in GLASD results.');
end

% ==== STACK, AVERAGE, AND STANDARD ERROR (GLASD ONLY) ==== 
stacked_glasd     = cat(3, valid_results{:});                     % 3D array 
mean_glasd        = mean(stacked_glasd, 3);                       % element-wise mean 
stderr_glasd      = std(stacked_glasd, 0, 3) / sqrt(size(stacked_glasd, 3));  % standard error

% --------------------------------------------------------------
% ==== READ & AGGREGATE NAIVE RESULTS OVER SAME VALID REPS ====
% --------------------------------------------------------------
num_valid   = sum(valid_idx);
naive_vals  = nan(num_valid, 2);  % columns: [RMSE, MAD]
cnt         = 0;

for DataRep = 1:NumDatasetReps
    if ~valid_idx(DataRep)
        continue;  % skip if GLASD file missing for this rep
    end
    
    cnt = cnt + 1;
    naive_filename = sprintf('Naive_rmse_mad_p_%d_n_%d_C_%s_Dist_%s_DataRep_%d.csv', ...
        p, n, Ctype, Dist, DataRep);
    
    if isfile(naive_filename)
        Nmat = readmatrix(naive_filename);
        naive_vals(cnt, :) = Nmat(1, end-1:end);   % last 2 cols = [RMSE, MAD]
    else
        warning('Naive file not found: %s', naive_filename);
        naive_vals(cnt, :) = [NaN NaN];
    end
end

valid_naive = ~any(isnan(naive_vals), 2);
naive_vals  = naive_vals(valid_naive, :);
if isempty(naive_vals)
    warning('No valid Naive results found; summary will exclude Naive.');
    include_naive = false;
else
    include_naive = true;
    mean_naive   = mean(naive_vals, 1);                          % 1 x 2
    stderr_naive = std(naive_vals, 0, 1) / sqrt(size(naive_vals, 1));
end

% --------------------------------------------------------------
% ==== READ & AGGREGATE TYLER RESULTS OVER SAME VALID REPS ====
% --------------------------------------------------------------
tyler_vals = nan(num_valid, 2);  % columns: [RMSE, MAD]
cnt = 0;

for DataRep = 1:NumDatasetReps
    if ~valid_idx(DataRep)
        continue;  % skip if GLASD file missing for this rep
    end
    
    cnt = cnt + 1;
    tyler_filename = sprintf('Tyler_rmse_mad_p_%d_n_%d_C_%s_Dist_%s_DataRep_%d.csv', ...
        p, n, Ctype, Dist, DataRep);
    
    if isfile(tyler_filename)
        Tmat = readmatrix(tyler_filename);
        tyler_vals(cnt, :) = Tmat(1, end-1:end);   % last 2 cols = [RMSE, MAD]
    else
        warning('Tyler file not found: %s', tyler_filename);
        tyler_vals(cnt, :) = [NaN NaN];
    end
end

valid_tyler = ~any(isnan(tyler_vals), 2);
tyler_vals  = tyler_vals(valid_tyler, :);
if isempty(tyler_vals)
    warning('No valid Tyler results found; summary will exclude Tyler.');
    include_tyler = false;
else
    include_tyler = true;
    mean_tyler   = mean(tyler_vals, 1);                          % 1 x 2
    stderr_tyler = std(tyler_vals, 0, 1) / sqrt(size(tyler_vals, 1));
end

% --------------------------------------------------------------
% ==== COMBINE: NAIVE, TYLER, THEN GLASD METHODS ====
% mean_glasd / stderr_glasd assumed [*, RMSE, MAD] (3 columns)
% --------------------------------------------------------------
rows_all   = {};
mean_all   = [];
stderr_all = [];

if include_naive
    rows_all   = [rows_all; {'Naive'}];
    mean_all   = [mean_all;   [NaN, mean_naive]];
    stderr_all = [stderr_all; [NaN, stderr_naive]];
end

if include_tyler
    rows_all   = [rows_all; {'Tyler'}];
    mean_all   = [mean_all;   [NaN, mean_tyler]];
    stderr_all = [stderr_all; [NaN, stderr_tyler]];
end

% GLASD-based methods, in fixed order:
glasd_method_names = {
    'Gaussian';
    'Huber';
    'Truncated';
    'Truncated(soft)';
    'Tukey'
};

rows_all   = [rows_all; glasd_method_names];
mean_all   = [mean_all;   mean_glasd];
stderr_all = [stderr_all; stderr_glasd];

% ==== ROUND NOW (AFTER COMBINING) ====
mean_all   = round(mean_all,   3);
stderr_all = round(stderr_all, 4);

% ==== FORMAT AS "mean (SE)" STRINGS ====
formatted_output = strings(size(mean_all, 1), 2);  % RMSE & MAD

for i = 1:size(mean_all, 1)
    for j = 2:3
        formatted_output(i, j-1) = sprintf('%.3f (%.4f)', ...
            mean_all(i, j), stderr_all(i, j));
    end
end

% ==== ADD TABLE HEADERS ====
summary_table = array2table(formatted_output, ...
    'VariableNames', {'RMSE', 'MAD'}, ...
    'RowNames', rows_all);

% ==== CONSTRUCT OUTPUT FILENAME ====
summary_filename = sprintf('Summary_rmse_mad_p_%d_n_%d_NumOptmzReps_%d_C_%s_Dist_%s.csv', ...
    p, n, NumOptmzReps, Ctype, Dist);

% ==== SAVE TO CSV WITH ROW NAMES ====
writetable(summary_table, summary_filename, 'WriteRowNames', true);
fprintf('Saved summary file: %s\n', summary_filename);
