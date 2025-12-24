clear; clc;

% ===================== USER INPUT =====================
files = {
    'C_naive.csv', ...
    'C_tyler.csv', ...
    'C_POET.csv', ...
    'C_huber_optimal.csv'
};

method_names = {
    'Naive', ...
    'Tyler', ...
    'POET', ...
    'GLASD-Huber'
};

thr = 0.3;   % edge threshold for density
% ======================================================

K = numel(files);

kappa_vals   = zeros(K,1);
lam_min_vals = zeros(K,1);
lam_max_vals = zeros(K,1);   % <<< NEW
mean_abs     = zeros(K,1);
median_abs   = zeros(K,1);
edge_density = zeros(K,1);

for k = 1:K

    % ---- Read matrix safely ----
    R = readmatrix(files{k});

    % ---- STRONG SAFETY CHECK ----
    if size(R,1) ~= size(R,2)
        error('Matrix %s is NOT square: %d x %d', ...
              files{k}, size(R,1), size(R,2));
    end

    % ---- Enforce symmetry (numerically essential) ----
    R = (R + R') / 2;

    % ---- Enforce unit diagonal SAFELY ----
    p = size(R,1);
    R(1:p+1:end) = 1;

    % ---- Eigenvalues (REAL enforce) ----
    eigvals = real(eig(R));

    lam_max = max(eigvals);
    lam_min = min(eigvals);

    lam_max_vals(k) = lam_max;   % <<< NEW STORAGE
    lam_min_vals(k) = lam_min;

    % ---- (1) CONDITION NUMBER (PD-aware) ----
    if lam_min <= 0
        kappa_vals(k) = Inf;
    else
        kappa_vals(k) = lam_max / lam_min;
    end

    % ---- Off-diagonal entries only ----
    off_diag = R(triu(true(p),1));

    % ---- (3) Mean |r_jk| ----
    mean_abs(k) = mean(abs(off_diag));

    % ---- (4) Median |r_jk| ----
    median_abs(k) = median(abs(off_diag));

    % ---- (5) Edge density ----
    edge_density(k) = mean(abs(off_diag) > thr);

end

% ===================== FINAL TABLE =====================
Results = table( ...
    kappa_vals, ...
    lam_min_vals, ...
    lam_max_vals, ...      % <<< ADDED
    mean_abs, ...
    median_abs, ...
    edge_density, ...
    'RowNames', method_names, ...
    'VariableNames', ...
    { ...
      'CondNumber', ...
      'Lambda_min', ...
      'Lambda_max', ...    % <<< ADDED
      'MeanAbsCorr', ...
      'MedianAbsCorr', ...
      'EdgeDensity_gt_0p3' ...
    } );

disp('===== REAL DATA COMPARISON METRICS =====');
disp(Results)

writetable(Results, ...
    'CRC_RobustCorr_Comparison_Table.csv', ...
    'WriteRowNames', true);

fprintf('\nSaved: CRC_RobustCorr_Comparison_Table.csv\n');



%% =========================================================
%   ORDERED EIGENVALUE SPECTRAL PLOT (LOG10 SCALE)
%   Methods: Naive, Tyler, GLASD-Huber (SKIP POET)
% ==========================================================

clearvars -except files method_names

use_idx     = [1 2 4];                      % Skip POET
plot_labels = method_names(use_idx);

figure; hold on;

for ii = 1:length(use_idx)

    k = use_idx(ii);

    % ---- Read correlation matrix ----
    R = readmatrix(files{k});

    % ---- Enforce symmetry and unit diagonal ----
    R = (R + R') / 2;
    p = size(R,1);
    R(1:p+1:end) = 1;

    % ---- Eigenvalues ----
    eigvals = sort(real(eig(R)), 'descend');

    % ---- Numerical safety ----
    eigvals(eigvals <= eps) = eps;

    % ---- LOG10 TRANSFORM ----
    logeig = log10(eigvals);

    % ---- Plot ----
    plot(logeig, 'LineWidth', 2);
end

legend(plot_labels, 'Location', 'southwest');
xlabel('Ordered Index');
ylabel('log_{10}(Eigenvalue)');
title('Ordered Eigenvalue Spectra (Log_{10} Scale)');
grid on;
box on;

set(gca, 'FontSize', 12);

saveas(gcf, 'Ordered_Eigenvalue_Log10_Spectral_Comparison.png');
fprintf('Saved: Ordered_Eigenvalue_Log10_Spectral_Comparison.png\n');
