function sumRes = L1_out_Cross_Validation(Model_Name, Avox, startx, bvals, qhat, starting_points)
%
% This function implements the Leave-One-Out Cross-Validation (LOO-CV) technique 
% to evaluate the generalization error of diffusion models such as the Ball-and-Stick 
% and Zeppelin-and-Stick models.
%
% INPUTS:
%   - Model_Name : String, name of the model function to be evaluated.
%   - Avox : Measured signal intensities for the voxel.
%   - startx : Initial parameter estimates for the optimization.
%   - bvals : Vector of b-values representing diffusion weightings.
%   - qhat : Matrix of gradient directions (Nx3).
%   - starting_points : Number of different starting points for the optimization.
%
% OUTPUT:
%   - sumRes : Vector of residual errors from each leave-one-out iteration

% Create the name of the function we are going to use to compute the value
% of the signal
signal_function = strcat(Model_Name, '_signal');

% Initialize
sumRes = inf(length(Avox), 1);

% Run the Monte-Carlo Cross-Validation algorithm
parfor i = 1 : length(Avox)

    % Select the index of the voxel we are going to use to train the model
    training_index = true(length(Avox), 1);
    training_index(i) = false;


    % Run the fitting procedure from different starting points
    if(min(Avox) > 0) % Accept only significant values

        if length(startx) ==5  % Ball-and-Stick and Zeppelin-and-Stick with Tortuosity models

            [best_parameters, ~, ~] = find_optimal_param_newdataset(Model_Name, startx, starting_points, Avox(training_index), bvals(1, training_index), qhat(:, training_index));

        else % Zeppelin-and-Stick model

            [best_parameters, ~, ~] = optimal_fit_ZeppelinStick(startx, starting_points, Avox(training_index), bvals(1, training_index), qhat(:, training_index));

        end
        
        % Calculate the signal with the parameters obtained from the
        % optimization process
        test_signal = feval(signal_function, best_parameters, bvals(1, ~training_index), qhat(:, ~training_index));

        % Compute the sum
        sumRes(i) = sum((Avox(~training_index) - test_signal').^2);
    end
end