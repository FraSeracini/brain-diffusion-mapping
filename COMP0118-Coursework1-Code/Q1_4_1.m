%% Q1.4.1

% Load and arrange the new data set
[D, bvals, qhat, TE] = Load_New_Dataset();

% Select a voxel
Avox = D(:, 1);

% Use the model of the Linear Diffusion Tensor to find a starting point for
% the optimization process 
startx = DT_starting_point(Avox, bvals, qhat);
num_runs = 15;

% Run the fitting procedure
[best_parameters, ~, ~] = find_optimal_param_newdataset('BallStickSSD_constraints', startx, num_runs, Avox, bvals, qhat);

% Evaluate the FIM at the maximum likelihood estimate of the Ball-and-Stick
% parameters
FIM = Fisher_Information_Matrix(best_parameters, Avox, bvals, qhat)

% Initialize the arrays where we are going to store the values for the A
% optimality and the T optimality for each shell
A_opt_value = zeros(36, 1);
T_opt_value = zeros(36, 1);

% Find the different shells 
shells = Find_Shells(bvals, TE);

% Calculate the A-optimality and T-optimality value for all the shells 
for i = 1 : 36
    A_opt_value(i) = A_Optimality(best_parameters, Avox(shells(i).indices), bvals(shells(i).indices), qhat(shells(i).indices));
    T_opt_value(i) = T_Optimality(best_parameters, Avox(shells(i).indices), bvals(shells(i).indices), qhat(shells(i).indices));
end

% Find the minimum A-optimality value
min_A = min(A_opt_value);

% Find the maximum T-optimality value
max_T = max(T_opt_value);

% Define the acceptable deviation from the minimum / maximum value to still consider it as a valid minimum / maximum.
tolerance = 0.000001;

% Select all the shells containing the minimum A-optimality value
valid_A_opt_shells = find(abs(A_opt_value - min_A) <= tolerance)

% Select all the shells containing the maximum T-optimality value
valid_T_opt_shells = find(abs(T_opt_value - max_T) <= tolerance)