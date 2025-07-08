%% Wild bootstrap

% Load and arrange the data
[dwis, qhat, bvals] = Load_DataSet_1();

% Choose the number of iterations
T = 1000;

% Select a voxel
Avox = dwis(:, 92, 65, 72);

% Impose the number of runs
num_runs = 15;

% Starting point
startx = [3.5e+00 3e-03 2.5e-01 0 0];

if(min(Avox) > 0) % Accept only significant values

    % Find the parameters associated to global minima of the objective function
    [best_parameters, min_RESNORM, ~] = find_optimal_parameters('BallStickSSD_constraints', startx, num_runs, Avox, bvals, qhat);

end

% Extract the estimated signal with the parameters obtained
A_est = BallStickSSD_constraints_signal(best_parameters, bvals, qhat);

% Calculate the residuals
residuals = Avox - A_est';

% Initialize the matrix
residual_samples = zeros(length(best_parameters), T);

% Apply Wild Bootstrap algorithm
for t = 1 : T

    % Select randomly the residuals
    resampled_residuals = datasample(residuals, 1);

    % Add the residuals to the model signal
    A_resampled = (A_est + resampled_residuals * rand)';

    if(min(A_resampled) > 0) % Accept only significant values

        % Find the parameters associated to global minima of the objective function
        [best_parameters, min_RESNORM, ~] = find_optimal_parameters('BallStickSSD_constraints', startx, num_runs, A_resampled, bvals, qhat);
    
    end

    % Store the values
    residual_samples(:, t) = best_parameters;
end

% Estimate the 2-sigma range and the 95% range for the parameters S0, d and f
sigma_ranges_wb = Two_Sigma_Range(residual_samples);
ninety_five_range_wb = NinetyFive_Range(residual_samples);
