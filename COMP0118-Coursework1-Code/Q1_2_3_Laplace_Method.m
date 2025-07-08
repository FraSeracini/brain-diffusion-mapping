%% Laplace's Method

% Load and arrange the data
[dwis, qhat, bvals] = Load_DataSet_1();

% Put data for one voxel into A
Avox = dwis(:, 92, 65, 72);

% Starting point
startx = [3.5e+00 3e-03 2.5e-01 0 0];

% Initialize different values
min_RESNORM = inf;
best_parameters = zeros(5);
final_Hessian = zeros(5, 5);
num_runs = 15;

% Run the fitting procedure from different starting points
for n = 1 : num_runs

    % Add normally distributed random numbers to obtain new starting points
    startx_n = add_randn_numbers(startx);
   
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
        [parameter_hat, RESNORM, ~, ~, ~, Hessian] = fminunc('BallStickSSD_constraints', startx_n, h, Avox, bvals, qhat);
    catch
        continue;
    end

    % Transformation to get the model parameters
    parameter_hat(1) = parameter_hat(1)^2;
    parameter_hat(2) = parameter_hat(2)^2;
    parameter_hat(3) = exp(-parameter_hat(3)^2);

    % Update the values if we find a RESNORM lower than the ones stored in min_RESNORM
    if(RESNORM < min_RESNORM)
        min_RESNORM = RESNORM;
        best_parameters = parameter_hat;
        final_Hessian = Hessian;
    end
end

% Compute sigma
sigma2 = (1 / (length(Avox) - length(best_parameters))) * min_RESNORM;

% Compute the covariance matrix
covariance = -2 * sigma2 * inv(-final_Hessian);

% Compute the standard deviation
std = sqrt(diag(covariance));

% Initialize the matrix where we are going to store the 2-sigma range
sigma_ranges = zeros(3, 2);

% Estimate the 2-sigma range for the parameters S0, d and f
for i = 1 : 3
    sigma_ranges(i, 1) = best_parameters(i) - 2 * std(i);
    sigma_ranges(i, 2) = best_parameters(i) + 2 * std(i);
end