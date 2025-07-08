%% Q1.1.4

% Repeat the fitting procedure in Q1_1_3 lots of times from	different starting
% points. To get each starting point, add normally distributed random
% numbers to the one we used in Q1_1_3

% Load and arrange the data
[dwis, qhat, bvals] = Load_DataSet_1();

% Select a voxel
Avox = dwis(:, 92, 65, 72);

% Number of different starting points
starting_points = 1000;

% Define a starting point for the non-linear fit
startx = [3.5e+00   3e-03   2.5e-01 0 0];

tic

% Run the fitting procedure from different starting points
if(min(Avox) > 0) % Accept only significant values

    % Find the parameters associated to global minima of the objective function
    [best_parameters, min_RESNORM, RESNORM_per_start_point] = find_optimal_parameters('BallStickSSD_constraints', startx, starting_points, Avox, bvals, qhat);

end

toc


% Define the acceptable deviation from the minimum value to still consider it as a valid minimum.
accepted_deviation = 0.1;

% Calculate the percentage of times the minimum RESNORM value was found
% across all runs from different starting points. Also, determine the
% number of runs required to be 95% confident in finding the global minimum.
[percentage, num_runs] = min_resnorm_percentage(RESNORM_per_start_point, starting_points, accepted_deviation)
