function X = simulate_data_from_corr(C, n, dist_type, varargin)
% SIMULATE_DATA_FROM_CORR generates n x p data with correlation matrix C
% dist_type: 'gaussian', 'gaussian_asymmetric_row', 'gaussian_asymmetric_col', 'gaussian_asymmetric_random', 't'

% --- Default values ---
default.RowContaminationProb = 0.1;
default.ColContaminationProb = 0.1;
default.RowShiftVal = 10;
default.ColShiftVal = 10;
default.RandomContaminationProb = 0.05;
default.RandomShiftVal = 100;
default.df = 3; % for t-dist

% --- Parse inputs ---
params = default;
for i = 1:2:length(varargin)
    params.(varargin{i}) = varargin{i+1};
end

% --- Core setup ---
p = size(C, 1);
L = chol(C, 'lower');

switch lower(dist_type)
    case 'gaussian'
        Z = randn(n, p);
        
    case 'gaussian_asymmetric_random'
        Z = randn(n, p);
        
        % Total number of elements to contaminate
        total_elements = n * p;
        num_contam = round(params.RandomContaminationProb * total_elements);
        
        % Randomly pick linear indices
        contam_inds = randperm(total_elements, num_contam);
        
        % Apply contamination
        Z(contam_inds) = Z(contam_inds) + params.RandomShiftVal;
        
    case 'gaussian_asymmetric_row'
        Z = randn(n, p);
        
        % Total number of rows to contaminate
        num_contam_rows = round(params.RowContaminationProb * n);
        
        % Randomly pick rows to contaminate
        contam_rows = randperm(n, num_contam_rows);
        
        for i = 1:num_contam_rows
            row_idx = contam_rows(i);
            
            % Randomly choose 30%–70% of columns to contaminate in this row
            prop = 0.3 + (0.7 - 0.3)*rand();
            num_cols = round(prop * p);
            
            % Randomly pick columns
            contam_cols = randperm(p, num_cols);
            
            % Apply contamination
            Z(row_idx, contam_cols) = Z(row_idx, contam_cols) + params.RowShiftVal;
        end
        
    case 'gaussian_asymmetric_col'
        Z = randn(n, p);
        
        % Total number of columns to contaminate
        num_contam_cols = round(params.ColContaminationProb * p);
        
        % Randomly pick columns to contaminate
        contam_cols = randperm(p, num_contam_cols);
        
        for i = 1:num_contam_cols
            col_idx = contam_cols(i);
            
            % Randomly choose 30%–70% of rows to contaminate in this column
            prop = 0.3 + (0.7 - 0.3) * rand();
            num_rows = round(prop * n);
            
            % Randomly pick rows
            contam_rows = randperm(n, num_rows);
            
            % Apply contamination
            Z(contam_rows, col_idx) = Z(contam_rows, col_idx) + params.ColShiftVal;
        end
        
    case 't'
        Z = randn(n, p);
        chi2_vals = chi2rnd(params.df, n, 1);
        Z = Z ./ sqrt(chi2_vals / params.df);
        
    otherwise
        error('Unknown distribution type: %s', dist_type);
end

X = Z * L';
end
