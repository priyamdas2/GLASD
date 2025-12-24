clear; clc;
addpath('./Output/');

%% ================= USER INPUTS =================
Methods = {'Naive','POET','FWZ','FWZ_latentF','GLASD'};

Regime  = 'random'; % 'diagonal' / 'random'
Setting = '2_componentwise';      % '1_elliptical' or '2_componentwise'
p       = 200;
nu      = 3;

dataDir = 'Output';
numReps = 10;
%% ===============================================

FinalMean = [];    % 9 × 4
FinalSE   = [];    % 9 × 4
MethodLabels = {};

for mm = 1:length(Methods)

    Method = Methods{mm};
    all_vals = [];

    for rep = 1:numReps

        filename = sprintf( ...
            '%s_errors_%s_setting%s_p_%d_nu_%d_rep_%d.csv', ...
            Method, Regime, Setting, p, nu, rep);

        fullpath = fullfile(dataDir, filename);

        if ~isfile(fullpath)
            error('Missing file: %s', filename);
        end

        X = readmatrix(fullpath);     % [1×K] or [M×K]
        X = X(:,1:4);                 % ---- KEEP ONLY FIRST 4 METRICS ----

        if rep == 1
            [numSubMethods, numMetrics] = size(X);
            all_vals = zeros(numSubMethods, numMetrics, numReps);
        else
            if ~isequal(size(X), [numSubMethods, numMetrics])
                error('Dimension mismatch in %s', filename);
            end
        end

        all_vals(:,:,rep) = X;
    end

    %% ---- Mean & Standard Error over reps ----
    MeanVals = mean(all_vals, 3);                       % [M × 4]
    SEVals   = std(all_vals, 0, 3) / sqrt(numReps);     % [M × 4]

    %% ---- Stack into final matrices ----
    FinalMean = [FinalMean; MeanVals];
    FinalSE   = [FinalSE;   SEVals];

    %% ---- Row labels ----
    if strcmp(Method,'GLASD')
        for k = 1:size(MeanVals,1)
            MethodLabels{end+1} = sprintf('GLASD_%d', k);
        end
    else
        MethodLabels{end+1} = Method;
    end
end

%% --------- CREATE TABLES ----------
MetricNames = {'Metric1','Metric2','Metric3','Metric4'};

MeanTable = array2table(FinalMean, ...
    'VariableNames', MetricNames, ...
    'RowNames', MethodLabels);

SETable = array2table(FinalSE, ...
    'VariableNames', MetricNames, ...
    'RowNames', MethodLabels);

disp('===== MEAN TABLE =====');
disp(MeanTable)

disp('===== STANDARD ERROR (SE) TABLE =====');
disp(SETable)

%% --------- SAVE TO CSV WITH SUMMARY NAMES ----------
MeanFile = sprintf('Output/Summary_MEAN_%s_%s_p%d_nu%d.csv', ...
                    Regime, Setting, p, nu);

SEFile   = sprintf('Output/Summary_SE_%s_%s_p%d_nu%d.csv', ...
                    Regime, Setting, p, nu);

writetable(MeanTable, MeanFile, 'WriteRowNames', true);
writetable(SETable,   SEFile,   'WriteRowNames', true);

fprintf('\nSaved files:\n%s\n%s\n', MeanFile, SEFile);


%% --------- NEW: COMBINED mean (SE) TABLE & FILE ----------
% Build a cell array of "mean (SE)" strings
[numRows, numCols] = size(FinalMean);
Formatted = cell(numRows, numCols);

for i = 1:numRows
    for j = 1:numCols
        m = FinalMean(i,j);
        s = FinalSE(i,j);

        if abs(m) >= 1e6
            % Large values: scientific notation
            Formatted{i,j} = sprintf('%.1e (%.1e)', m, s);
        else
            % Moderate values: fixed decimals
            Formatted{i,j} = sprintf('%.1f (%.2f)', m, s);
        end
    end
end

MeanSETable = cell2table(Formatted, ...
    'VariableNames', MetricNames, ...
    'RowNames', MethodLabels);

MeanSEFile = sprintf('Output/Summary_MeanSE_%s_%s_p%d_nu%d.csv', ...
                     Regime, Setting, p, nu);

writetable(MeanSETable, MeanSEFile, 'WriteRowNames', true);

fprintf('%s\n', MeanSEFile);