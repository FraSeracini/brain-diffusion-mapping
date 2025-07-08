function [best_parameters, min_RESNORM, RESNORM_per_start_point] = optimal_fit_Rician(startx, num_runs, Avox, bvals, qhat)
%
% This function performs multiple optimization runs to identify the global minimum
% of the objective function for the Ball-and-Stick model, considering the effects 
% of Rician noise. It returns the optimal parameters, the minimum residual sum of 
% squares (RESNORM), and an array of RESNORM values from each optimization run.
%
% INPUTS:
%   - startx : Initial starting point for the optimization.
%   - num_runs : Number of different starting points for multiple optimizations.
%   - Avox : Measured signal intensities for the voxel.
%   - bvals : Vector of b-values representing diffusion weightings.
%   - qhat : Matrix of gradient directions (Nx3).
%
% OUTPUTS:
%   - best_parameters : Parameters associated with the global minimum.
%   - min_RESNORM : Minimum residual sum of squares (RESNORM) obtained.
%   - RESNORM_per_start_point : Array containing RESNORM values from each optimization run.

% Initialize the array where we are going to store all the RESNORM values
RESNORM_per_start_point  = zeros(1, num_runs);

% Initialize the minimum RESNORM and the set of model paramaeters assocciated
min_RESNORM = inf;
best_parameters = zeros(6);

% Run the fitting procedure from different starting points
for n = 1 : num_runs

    % Add normally distributed random numbers to obtain new starting points
    startx_n(1) = startx(1) + randn*1e+03;
    startx_n(2) = startx(2) + randn*1e-03;
    startx_n(3) = startx(3) + randn*1e-01;
    startx_n(4) = startx(4) + randn*1e-01;
    startx_n(5) = startx(5) + randn*1e-01;

    % Impose that S0 >= 0
    while (startx_n(1) < 0)
        startx_n(1) = startx(1) + randn*1e+03;
    end

    % Impose that d >= 0
    while (startx_n(2) < 0)
        startx_n(2) = startx(2) + randn*1e-03;
    end

    % Impose that 0 <= f <= 1
    while startx_n(3) < 0 || startx_n(3) >1
        startx_n(3) = startx(3) + randn*1e-01;
    end

    % Inverse transformation to mantain the same starting point
    startx_n(1) = sqrt(startx_n(1));
    startx_n(2) = sqrt(startx_n(2));
    startx_n(3) = sqrt(-log(startx_n(3)));

    % Define various options for the non-linear fitting algorithm
    h=optimset('MaxFunEvals',20000,...
        'Algorithm','quasi-newton',...
        'MaxIter',2000,...
        'TolX',1e-10,...
        'TolFun',1e-10,...
        'Display','off');

    % Now run the fitting
    try
        [parameter_hat, RESNORM, ~, ~] = fminunc('BallStickSSD_Rician', startx_n, h, Avox, bvals, qhat);
    catch
        continue;
    end

    % Transformation to get the model parameters
    parameter_hat(1) = parameter_hat(1)^2;
    parameter_hat(2) = parameter_hat(2)^2;
    parameter_hat(3) = exp(-parameter_hat(3)^2);

    % Store the RESNORM value in the array
    RESNORM_per_start_point(1, n) = RESNORM;

    % Update the values if we find a RESNORM lower than the
    % ones stored in min_RESNORM
    if(RESNORM < min_RESNORM)
        min_RESNORM = RESNORM;
        best_parameters = parameter_hat;
    end
end
end