%% Q1.3.1

% Adapt	the	code from Q1.1.3 to	fit	the	ball and stick model to	this new
% dataset

% Load and arrange the new data set
[D, bvals, qhat] = Load_New_Dataset();

% Number of different starting points
starting_points = 1000;

% Initialize the array where we are going to store all the RESNORM values
RESNORM_per_start_point  = zeros(1, starting_points);

% Initialize the minimum RESNORM and the set of model paramaeters assocciated
min_RESNORM_BS = inf;
best_parameters_BS = zeros(5);

% Select a voxel
Avox = D(:, 1);


% Use the model of the Linear Diffusion Tensor to find a starting point for
% the optimization process 
startx = DT_starting_point(Avox, bvals, qhat);

tic
% Run the fitting procedure from different starting points
if(min(Avox) > 0) % Accept only significant values

    % Find the parameters associated to global minima of the objective function
    [best_parameters_BS, min_RESNORM_BS, RESNORM_per_start_point] = find_optimal_param_newdataset('BallStickSSD_constraints', startx, starting_points, Avox, bvals, qhat);

end
toc

best_parameters_BS
min_RESNORM_BS

% Define the acceptable deviation from the minimum value to still consider it as a valid minimum.
accepted_deviation = 0.001;

% Calculate the percentage of times the minimum RESNORM value was found
% across all runs from different starting points. Also, determine the
% number of runs required to be 95% confident in finding the global minimum.
[percentage, num_runs] = min_resnorm_percentage( RESNORM_per_start_point, starting_points, accepted_deviation)

% Extract the estimated signal with the parameters obtained
Aest = BallStickSSD_constraints_signal(best_parameters_BS, bvals, qhat);

figure;
plot(Avox, 'bs')
hold on;
plot(Aest, 'rx')
xlabel('k');
ylabel('S');
legend('Data', 'Model');

% Save the results
save('min_RESNORM_BS','min_RESNORM_BS');
