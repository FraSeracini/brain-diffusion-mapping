%% Q1.1.5

% Create parameter maps	over one slice of the image volume. Produce maps of
% S0, d, f and the residual error RESNORM. Then produce a map of the
% fibre direction n

% Load and arrange the data
[dwis, qhat, bvals] = Load_DataSet_1();

% Define a starting point for the non-linear fit
startx = [3.5e+00   3e-03   2.5e-01 0 0];

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
num_runs = 15;

% Run the fitting procedure for each voxel in the slice
for i = 1 : 145
    for j = 1 : 174

        % Update the voxel in the chosen slice
        Avox = dwis(:, i, j, 72);

        if(min(Avox) > 0) % Accept only significant values

            % Find the parameters associated to global minima of the objective function
            [best_parameters, min_RESNORM, ~] = find_optimal_parameters('BallStickSSD_constraints', startx, num_runs, Avox, bvals, qhat);

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