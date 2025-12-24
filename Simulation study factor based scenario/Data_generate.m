%% ============================================
% Fan–Wang–Zhong (2019) Data Generator
% BEST COMPARISON DESIGN:
%   - Three Sigma_u regimes:
%       'sparse'   -> banded AR(1) 
%       'diagonal' -> pure diagonal 
%       'random'   -> random sparse graph 
%   - Tunable diagonal scale sigma_u_diag
% Output filenames encode regime to avoid confusion.
% ============================================

clear; clc;
addpath('./FanEtAl2019/');
addpath('./GLASD/');
addpath('./supp funs/');

% ------------------- Basic settings -------------------
p   = 200;
r   = 3;
n   = p / 2;
nu  = 3;
num_exp = 10;
for rep = 1:num_exp     % Monte Carlo replicate index
    
    rng(rep);
    
    % Noise diagonal scale (std of idiosyncratic components)
    sigma_u_diag = 2;
    
    % Sparse structure parameters
    rho_off   = 0.3;   % band strength for AR(1)
    edge_prob = 0.05;  % density for random sparse graph (~5%)
    
    % We now generate data for THREE regimes
    Sigma_u_regime_list = {'sparse', 'diagonal', 'random'};
    
    % --------------------------------------------
    % Generate B once (fixed across all datasets)
    % --------------------------------------------
    B = randn(p, r);
    fname_B = sprintf('Data/B_matrix_p_%d_nu_%d_rep_%d.csv', p, nu, rep);
    writematrix(B, fname_B);
    
    % ======================================================
    % Loop over Sigma_u regimes
    % ======================================================
    for g = 1:numel(Sigma_u_regime_list)
        
        Sigma_u_regime = Sigma_u_regime_list{g};
        
        % --------------------------------------------
        % Construct Sigma_u for this regime
        % --------------------------------------------
        switch lower(Sigma_u_regime)
            
            case 'diagonal'
                % GLASD-favorable: pure diagonal
                Sigma_u_true = (sigma_u_diag^2) * eye(p);
                
            case 'sparse'
                % FWZ-favorable: diagonally dominant banded AR(1)
                Sigma_u_true = (sigma_u_diag^2) * eye(p);
                for i = 1:p-1
                    Sigma_u_true(i,   i+1) = rho_off;
                    Sigma_u_true(i+1, i  ) = rho_off;
                end
                
            case 'random'
                % General random sparse graph
                Sigma_u_true = zeros(p);
                
                for i = 1:p-1
                    for j = i+1:p
                        if rand < edge_prob
                            val = rho_off * sign(randn);
                            Sigma_u_true(i,j) = val;
                            Sigma_u_true(j,i) = val;
                        end
                    end
                end
                
                % Enforce diagonal dominance for SPD
                for i = 1:p
                    Sigma_u_true(i,i) = sigma_u_diag^2 + sum(abs(Sigma_u_true(i,:)));
                end
                
            otherwise
                error('Unknown Sigma_u_regime: %s.', Sigma_u_regime);
        end
        
        % Save Sigma_u_true
        fname_Su = sprintf('Data/SigmaU_%s_p_%d_nu_%d_rep_%d.csv', ...
            Sigma_u_regime, p, nu, rep);
        writematrix(Sigma_u_true, fname_Su);
        
        % Cholesky for noise generation
        Lu = chol(Sigma_u_true, 'lower');
        
        %% ============================================
        % Setting 1: Elliptical multivariate t
        % ============================================
        
        % Factors: elliptical t with identity covariance
        F1 = trnd(nu, n, r);   % n x r
        
        % Idiosyncratic noise
        E1 = trnd(nu, n, p);
        U1 = E1 * Lu';
        
        % Observations
        Y1 = F1 * B' + U1;
        
        fname_Y1 = sprintf('Data/Y_%s_setting1_elliptical_p_%d_nu_%d_rep_%d.csv', ...
            Sigma_u_regime, p, nu, rep);
        fname_F1 = sprintf('Data/F_%s_setting1_elliptical_p_%d_nu_%d_rep_%d.csv', ...
            Sigma_u_regime, p, nu, rep);
        
        writematrix(Y1, fname_Y1);
        writematrix(F1, fname_F1);
        
        %% ============================================
        % Setting 2: Componentwise iid t
        % ============================================
        
        F2 = trnd(nu, n, r);
        
        E2 = trnd(nu, n, p);
        U2 = E2 * Lu';
        
        Y2 = F2 * B' + U2;
        
        fname_Y2 = sprintf('Data/Y_%s_setting2_componentwise_p_%d_nu_%d_rep_%d.csv', ...
            Sigma_u_regime, p, nu, rep);
        fname_F2 = sprintf('Data/F_%s_setting2_componentwise_p_%d_nu_%d_rep_%d.csv', ...
            Sigma_u_regime, p, nu, rep);
        
        writematrix(Y2, fname_Y2);
        writematrix(F2, fname_F2);
        
        fprintf('Generated data for Sigma_u_regime = %s\n', Sigma_u_regime);
    end
    
end
