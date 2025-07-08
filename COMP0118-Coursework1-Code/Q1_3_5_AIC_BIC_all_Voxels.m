%% Q1.3.5
% Check if the results obtained in Q1.3.3 are consistents across the other
% five voxels

% Load and arrange the new data set
[D, bvals, qhat] = Load_New_Dataset();

% Initialize the values
SSD_DT = zeros(1, 5);
min_RESNORM_BS = inf(1, 5);
min_RESNORM_ZS = inf(1, 5);
min_RESNORM_ZST = inf(1, 5);
AIC_DT = zeros(1, 5);
AIC_BS = zeros(1, 5);
AIC_ZS = zeros(1, 5);
AIC_ZST = zeros(1, 5);
BIC_DT = zeros(1, 5);
BIC_BS = zeros(1, 5);
BIC_ZS = zeros(1, 5);
BIC_ZST = zeros(1, 5);

% Number of parameters in each model
n_parameters_DT = 7;
n_parameters_BS = 5;
n_parameters_ZS = 6;
n_parameters_ZST = 5;

% Apply a correction to avoid errors
D(1338, 5) = 0.1;

for i = 2 : 6

    % Select a voxel
    Avox = D(:, i);

    % Diffusion Tensor model
    [x, dt, S, ~] = Linear_Diffusion_Tensor(Avox, bvals, qhat);

    % Calculate the SSD for the Diffusion Tensor Model
    SSD_DT(1, i-1) = sum((Avox - S).^2);

    % Use the model of the Linear Diffusion Tensor to find a starting point for
    % the optimization process
    startx = DT_starting_point(Avox, bvals, qhat);
    startx_ZS = [startx startx(2)];

    % Number of different starting points
    starting_points = 100;

    % Initialize the array where we are going to store all the RESNORM values
    RESNORM_per_start_point  = zeros(1, starting_points);

    % Initialize the minimum RESNORM and the set of model paramaeters assocciated
    best_parameters_BS = zeros(5);

    % Run the fitting procedure from different starting points
    if(min(Avox) > 0) % Accept only significant values

        % Find the parameters associated to global minima of the objective function
        [best_parameters_BS, min_RESNORM_BS(1, i-1), RESNORM_per_start_point] = find_optimal_param_newdataset('BallStickSSD_constraints', startx, starting_points, Avox, bvals, qhat);

    end

    %Zeppelin-and-Stick model

    % Number of different starting points
    starting_points = 100;

    if min(Avox) > 0 % Accept only significant values

        % Find the parameters associated to global minima of the objective function
        [best_parameters_ZS, min_RESNORM_ZS(1, i-1), RESNORM_per_start_point] = optimal_fit_ZeppelinStick(startx_ZS, starting_points, Avox, bvals, qhat);

    end

    % Zeppelin-and-Stick with tortuosity model

    if(min(Avox) > 0) % Accept only significant values

        % Find the parameters associated to global minima of the objective function
        [best_parameters_ZST, min_RESNORM_ZST(1, i-1), RESNORM_per_start_point] = find_optimal_param_newdataset('ZeppelinStickTortuosity', startx, starting_points, Avox, bvals, qhat);

    end

    % Compute the AIC for each model
    AIC_DT(1, i-1) = Compute_AIC(n_parameters_DT, Avox, SSD_DT(1, i-1));
    AIC_BS(1, i-1) = Compute_AIC(n_parameters_BS, Avox, min_RESNORM_BS(1, i-1));
    AIC_ZS(1, i-1) = Compute_AIC(n_parameters_ZS, Avox, min_RESNORM_ZS(1, i-1));
    AIC_ZST(1, i-1) = Compute_AIC(n_parameters_ZST, Avox, min_RESNORM_ZST(1, i-1));

    % Compute the AIC for each model
    BIC_DT(1, i-1) = Compute_BIC(n_parameters_DT, Avox, SSD_DT(1, i-1));
    BIC_BS(1, i-1) = Compute_BIC(n_parameters_BS, Avox, min_RESNORM_BS(1, i-1));
    BIC_ZS(1, i-1) = Compute_BIC(n_parameters_ZS, Avox, min_RESNORM_ZS(1, i-1));
    BIC_ZST(1, i-1) = Compute_BIC(n_parameters_ZST, Avox, min_RESNORM_ZST(1, i-1));
end
