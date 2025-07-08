%% Q1.3.2

% Adapt the code	further	to fit various different models: Diffusion tensor
% model, Zeppelin-and-Stick model, Zeppelin-and-Stick with tortuosity model

% Load and arrange the new data set
[D, bvals, qhat] = Load_New_Dataset();

% Select a voxel
Avox = D(:, 1);

%% Diffusion Tensor model

[x, dt, S, ~] = Linear_Diffusion_Tensor(Avox, bvals, qhat);

% Calculate the SSD for the Diffusion Tensor Model
SSD_DT = sum((Avox - S).^2)

% Save value
save('min_RESNORM_DT','SSD_DT');

% Visualize the fit to the data by comparing the predictions from the model
% with the fitted parameters obtained to the actual measurements
figure;
plot(Avox, 'bs')
hold on;
plot(S, 'rx')
xlabel('k');
ylabel('S');
legend('Data', 'Model');

%% Zeppelin-and-Stick model

% Use the model of the Linear Diffusion Tensor to find a starting point for
% the optimization process 
startx = DT_starting_point(Avox, bvals, qhat);
startx_ZS = [startx startx(2)];

% Select the number of starting points
starting_points = 1000;

tic

if min(Avox) > 0 % Accept only significant values
    
    % Find the parameters associated to global minima of the objective function
    [best_parameters_ZS, min_RESNORM_ZS, RESNORM_per_start_point] = optimal_fit_ZeppelinStick(startx_ZS, starting_points, Avox, bvals, qhat);

end
toc

% Define the acceptable deviation from the minimum value to still consider it as a valid minimum.
accepted_deviation = 0.001;

% Calculate the percentage of times the minimum RESNORM value was found
% across all runs from different starting points. Also, determine the
% number of runs required to be 95% confident in finding the global minimum.
[percentage_ZS, num_runs_ZS] = min_resnorm_percentage(RESNORM_per_start_point, starting_points, accepted_deviation);

% Save value
save('min_RESNORM_ZS','min_RESNORM_ZS');

% Extract the estimated signal with the parameters obtained
Aest = ZeppelinStickSSD_signal(best_parameters_ZS, bvals, qhat);

% Visualize the fit to the data by comparing the predictions from the model
% with the fitted parameters obtained to the actual measurements
figure;
plot(Avox, 'bs')
hold on;
plot(Aest, 'rx')
xlabel('k');
ylabel('S');
legend('Data', 'Model');

%% Zeppelin-and-Stick with tortuosity model

tic
if(min(Avox) > 0) % Accept only significant values

    % Find the parameters associated to global minima of the objective function
    [best_parameters_ZST, min_RESNORM_ZST, RESNORM_per_start_point] = find_optimal_param_newdataset('ZeppelinStickTortuosity', startx, starting_points, Avox, bvals, qhat);

end
toc

% Define the acceptable deviation from the minimum value to still consider it as a valid minimum.
accepted_deviation = 0.001;

% Calculate the percentage of times the minimum RESNORM value was found
% across all runs from different starting points. Also, determine the
% number of runs required to be 95% confident in finding the global minimum.
[percentage_ZST, num_runs_ZST] = min_resnorm_percentage( RESNORM_per_start_point, starting_points, accepted_deviation)

% Save value
save('min_RESNORM_ZST','min_RESNORM_ZST');

% Extract the estimated signal with the parameters obtained
Aest = ZeppelinStickTortuosity_signal(best_parameters_ZST, bvals, qhat);

% Visualize the fit to the data by comparing the predictions from the model
% with the fitted parameters obtained to the actual measurements
figure;
plot(Avox, 'bs')
hold on;
plot(Aest, 'rx')
xlabel('k');
ylabel('S');
legend('Data', 'Model');