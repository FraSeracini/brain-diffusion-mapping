%% Q1.1.6

% Load and arrange the data
[dwis, qhat, bvals] = Load_DataSet_1();

%% fmincon

% An alternative way to	constrain your parameter estimates to realistic	ranges, as	
% opposed to the transformation	method in Q1.1.2, is to	use	constrained	
% optimization via the fmincon function.

% Select a voxel
Avox = dwis(:, 92, 65, 72);

% Define a starting point for the non-linear fit
startx = [3.5e+00   3e-03   2.5e-01 0 0];

% Number of different starting points
starting_points = 1000;

% Initialize the array where we are going to store all the RESNORM values
RESNORM_per_start_point  = zeros(1, starting_points);

% Initialize the minimum RESNORM and the set of model paramaeters assocciated
min_RESNORM = inf;
best_parameters = zeros(5);

% Constraints
lb = [0, 0, 0, -inf, -inf];
ub = [inf, inf, 1, inf, inf];

tic

% Run the fitting procedure from multiple starting points to understand the
% minimum number of runs we have to perform to be 95% sure to find the
% global minima
if min(Avox) > 0 % Accept only significant values
    
    % Run the fitting procedure from different starting points
    for n = 1 : starting_points

        % Add normally distributed random numbers to obtain new starting points
        startx_n = add_randn_numbers(startx);

        % Define various options for the non-linear fitting algorithm.
        h=optimoptions('fmincon',...
            'MaxFunctionEvaluations',20000,...
            'Algorithm','sqp',...
            'MaxIter', 2000,...
            'TolX',1e-10,...
            'TolFun',1e-10,...
            'Display','off');

        % Define the objective function
        objFun = @(x) BallStickSSD(x, Avox, bvals, qhat);

        % Now run the fitting
        try
            [parameter_hat, RESNORM, EXITFLAG, OUTPUT] = fmincon(objFun, startx_n, [], [], [], [], lb, ub, [], h);
        catch
            continue
        end

        % Store the RESNORM value in the array
        RESNORM_per_start_point(1, n) = RESNORM;

        % Update the values if we find a RESNORM lower than the ones stored in min_RESNORM
        if(RESNORM < min_RESNORM)
            min_RESNORM = RESNORM;
            best_parameters = parameter_hat;
        end

    end
end

toc

% Define the acceptable deviation from the minimum value to still consider it as a valid minimum.
accepted_deviation = 0.1;

% Calculate the percentage of times the minimum RESNORM value was found
% across all runs from different starting points. Also, determine the
% number of runs required to be 95% confident in finding the global minimum.
[percentage, num_runs] = min_resnorm_percentage(RESNORM_per_start_point, starting_points, accepted_deviation)

%% 

% Run the procedure over all the voxel in one slice to map the different
% parameteres

% Define a starting point for the non-linear fit
startx = best_parameters;

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

for i = 1 : 145
    for j = 1 : 174
        
        min_RESNORM = inf;
        best_parameters = zeros(1, 5);
       
        % Update the voxel in the chosen slice
        Avox=dwis(:, i, j, 72);

        if min(Avox) > 0 % Accept only significant values
            
            % Run the fitting procedure from different starting points
            for n = 1 : num_runs

                % Add normally distributed random numbers to obtain new starting points
                startx_n = add_randn_numbers(startx);

                % Define various options for the non-linear fitting algorithm.
                h=optimoptions('fmincon',...
                    'MaxFunctionEvaluations',20000,...
                    'Algorithm','sqp',...
                    'MaxIter', 2000,...
                    'TolX',1e-10,...
                    'TolFun',1e-10,...
                    'Display','off');

                % Define the objective function
                objFun = @(x) BallStickSSD(x, Avox, bvals, qhat);

                % Run the optimization
                try
                    [parameter_hat, RESNORM, EXITFLAG, OUTPUT] = fmincon(objFun, startx_n, [], [], [], [], lb, ub, [], h);
                catch
                    continue
                end
                
                
                % Update the values if we find a RESNORM lower than the
                % ones stored in min_RESNORM
                if(RESNORM < min_RESNORM)
                    min_RESNORM = RESNORM;
                    best_parameters = parameter_hat;
                end

            end

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

