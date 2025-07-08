%% Parametric Bootstrap

% Load and arrange the data
[dwis, qhat, bvals] = Load_DataSet_1();

% Starting point
startx = [3.5e+00 3e-03 2.5e-01 0 0];

% Select a voxel
Avox = dwis(:, 92, 65, 72);

% Impose the number of runs
num_runs = 15;

if(min(Avox) > 0) % Accept only significant values

    % Find the parameters associated to global minima of the objective function
    [best_parameters, min_RESNORM, ~] = find_optimal_parameters('BallStickSSD_constraints', startx, num_runs, Avox, bvals, qhat);

end

% Choose the number of iterations
T = 1000;

% Extract the estimated signal with the parameters obtained
signal = BallStickSSD_constraints_signal(best_parameters, bvals, qhat);

% Compute the standard deviation
sigma = sqrt((1 / (length(Avox) - length(best_parameters))) * sum(Avox - (signal)').^2);

% Initialize the matrix where we are going to store the samples
total_samples = zeros(length(best_parameters), T);

% Apply the parametric bootstrap algorithm
for t = 1 : T

    % Create samples from noise distribution N(0, sigma)
    E_samples = normrnd(0, sigma,length(Avox), 1);

    % Synthesize bootstrap data set
    A_hat = (signal)' + E_samples;

    % Inverse transformation to mantain the same starting point
    startx_n = startx;
    startx_n(1) = sqrt(startx(1));
    startx_n(2) = sqrt(startx(2));
    startx_n(3) = sqrt(-log(startx(3)));

    % Define various options for the non-linear fitting algorithm
    h=optimset('MaxFunEvals',500,...
        'Algorithm','quasi-newton',...
        'MaxIter',200,...
        'TolX',1e-10,...
        'TolFun',1e-10,...
        'Display','off');

    % Now run the fitting
    try
        [parameter_hat, RESNORM, ~, ~] = fminunc('BallStickSSD_constraints', startx_n, h, A_hat, bvals, qhat);
    catch
        continue;
    end

    % Transformation to get the model parameters
    parameter_hat(1) = parameter_hat(1)^2;
    parameter_hat(2) = parameter_hat(2)^2;
    parameter_hat(3) = exp(-parameter_hat(3)^2);

    % Store the values
    total_samples(:, t) = parameter_hat;
end

% Estimate the 2-sigma range and the 95% range for the parameters S0, d and f
sigma_ranges_pb = Two_Sigma_Range(total_samples);
ninety_five_range_pb = NinetyFive_Range(total_samples);
