function sumRes = MCarlo_Cross_Validation(Model_Name, Avox, test_dimension, num_runs, startx, bvals, qhat, starting_points)
%
% This function implements the Monte-Carlo Cross-Validation (MCCV) technique 
% to assess the predictive performance of diffusion models such as the Ball-and-Stick 
% and Zeppelin-and-Stick models.
%
% INPUTS:
%   - Model_Name : String, name of the model function to be evaluated.
%   - Avox : Measured signal intensities for the voxel.
%   - test_dimension : Number of randomly selected test samples per run.
%   - num_runs : Number of Monte-Carlo iterations.
%   - startx : Initial parameter estimates for the optimization.
%   - bvals : Vector of b-values representing diffusion weightings.
%   - qhat : Matrix of gradient directions (Nx3).
%   - starting_points : Number of different starting points for the optimization.
%
% OUTPUT:
%   - sumRes : Vector of residual errors from each Monte-Carlo iteration.

% Create the name of the function we are going to use to compute the value
% of the signal
signal_function = strcat(Model_Name, '_signal');

% Initiliaze
sumRes = inf(num_runs, 1);

% Run the Monte-Carlo Cross-Validation algorithm
for i = 1 : num_runs
    
    % Select the index of the voxel we are going to use to train the model
    training_index = true(length(Avox), 1);
    test_samples = randperm(length(Avox), test_dimension);
    training_index(test_samples) = false;

    if(min(Avox) > 0) % Accept only significant values

        % Run the fitting procedure
        if length(startx) ==5 % Ball-and-Stick and Zeppelin-and-Stick with Tortuosity models
            
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