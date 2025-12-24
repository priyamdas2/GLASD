% ==== SETUP ====
Ctypes = {'random', 'SparseUniform', 'BlockToeplitz'};
pn_pairs = [20 100; 50 500];
Dists = {
    'gaussian_asymmetric_row',    'Row Contamination';
    'gaussian_asymmetric_col',    'Column Contamination';
    'gaussian_asymmetric_random', 'Elementwise Contamination';
    't',                           't-distribution'
};
method_names = {'Gaussian', 'Huber', 'Truncated', 'Tukey'};

% ==== PREPARE STORAGE ====
row_labels = strings(24, 1);
RMSE_matrix = strings(24, 4);  % 24 rows × 4 methods

row_idx = 1;
for d = 1:size(Dists,1)
    dist_code = Dists{d,1};
    dist_label = Dists{d,2};

    for c = 1:length(Ctypes)
        for pni = 1:size(pn_pairs, 1)
            p = pn_pairs(pni, 1);
            n = pn_pairs(pni, 2);
            data_label = sprintf('C=%s, p=%d, n=%d (%s)', ...
                Ctypes{c}, p, n, dist_label);
            row_labels(row_idx) = data_label;

            % Read summary file
            fname = sprintf('Summary_rmse_mad_p_%d_n_%d_NumOptmzReps_10_C_%s_Dist_%s.csv', ...
                p, n, Ctypes{c}, dist_code);

            if isfile(fname)
                T = readtable(fname, 'ReadRowNames', true);
                for m = 1:length(method_names)
                    method = method_names{m};
                    if ismember(method, T.Properties.RowNames)
                        RMSE_matrix(row_idx, m) = string(T{method, 'RMSE'});
                    else
                        RMSE_matrix(row_idx, m) = "NA";
                    end
                end
            else
                warning('Missing file: %s', fname);
                RMSE_matrix(row_idx, :) = "NA";
            end

            row_idx = row_idx + 1;
        end
    end
end

% ==== CONVERT TO TABLE ====
final_table_24x4 = array2table(RMSE_matrix, ...
    'VariableNames', method_names, ...
    'RowNames', row_labels);

% ==== SAVE TO CSV ====
writetable(final_table_24x4, 'Final_RMSE_24x4_Summary.csv', 'WriteRowNames', true);
fprintf('Saved final 24x4 summary table to: Final_RMSE_24x4_Summary.csv\n');
