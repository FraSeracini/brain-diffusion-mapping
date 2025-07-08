%% Q1.1.7
% Alternative objective	function for the Rician	noise model

%% Parameter estimate

% Load and arrange the data
[dwis, qhat, bvals] = Load_DataSet_1();

% Select a voxel
Avox = dwis(:, 92, 65, 72);

% Define a starting point for the non-linear fit
startx = [3.5e+00   3e-03   2.5e-01 0 0 200];

% Number of starting points
starting_points = 1000;


% Run the fitting procedure from different starting points
if(min(Avox) > 0) % Accept only significant values
    
    % Find the parameters associated to global minima of the objective function
    [best_parameters, min_RESNORM, RESNORM_per_start_point] = optimal_fit_Rician(startx, starting_points, Avox, bvals, qhat);

end
best_parameters
min_RESNORM
% Define the acceptable deviation from the minimum value to still consider it as a valid minimum.
accepted_deviation = 0.1;

% Calculate the percentage of times the minimum RESNORM value was found
% across all runs from different starting points. Also, determine the
% number of runs required to be 95% confident in finding the global minimum
[percentage, num_runs] = min_resnorm_percentage(RESNORM_per_start_point, starting_points, accepted_deviation)

%% Parameters maps

% Create parameter maps	over one slice of the image volume. Produce maps of
% S0, d, f and the residual error RESNORM. Then produce a map of the
% fibre direction n

% Select a voxel
Avox = dwis(:, 92, 65, 72);

% Use the model of the Linear Diffusion Tensor to find a starting point for
% the optimization process 
startx = DT_starting_point(Avox, bvals, qhat);

tic
% Initialize the different maps
S0_map = zeros(145, 174);
d_map = zeros(145, 174);
f_map = zeros(145, 174);
theta_map = zeros(145, 174);
phi_map = zeros(145, 174);
fibre_direction_map = zeros(145, 174, 2);
RESNORM_map = zeros(145, 174);

% Impose the number of runs
num_runs = 10;

% Run the fitting procedure for each voxel in the slice
parfor i = 1 : 145
    for j = 1 : 174

        % Update the voxel in the chosen slice
        Avox=dwis(:, i, j, 72);

        if(min(Avox) > 0) % Accept only significant values

            % Find the parameters associated to global minima of the objective function
            [best_parameters, min_RESNORM, ~] = optimal_fit_Rician(startx, num_runs, Avox, bvals, qhat);

            % Store the values of the parameters associated with the global
            % minima in the different maps
            S0_map(i, j) = best_parameters(1);
            d_map(i, j) = best_parameters(2);
            f_map(i, j) = best_parameters(3);
            theta_map(i, j) = best_parameters(4);
            phi_map(i, j) = best_parameters(5);
            RESNORM_map(i, j) = min_RESNORM;
        end
    end
end

% Compute the fibre direction map
for i = 1 : 145
    for j = 1 : 174
        fibre_direction_map(i, j, :) = f_map(i, j) * [cos(phi_map(i, j)) * sin(theta_map(i, j)) sin(phi_map(i, j)) * sin(theta_map(i, j))];
    end
end

% Correction of out-of-scale values
RESNORM_map(RESNORM_map == inf) = 0;
toc

% Visualize the maps
figure
imshow(flipud(S0_map'), [])
title('S0 Map')

figure
imshow(flipud(d_map'), [])
title('d Map')

figure
imshow(flipud(f_map'), [])
title('f Map')

figure
imshow(flipud(RESNORM_map'), [])
title('RESNORM map')

fig=figure;
set(fig,'Position', [0, 0, 300, 350])
quiver(flipud(fibre_direction_map(:, :, 1)'), flipud(fibre_direction_map(:, :, 2)'))
title('Fibre Direction Map')
