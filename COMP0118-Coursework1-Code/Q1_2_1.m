%% Q1.2.1

% Use the classical	bootstrap procedure	to estimate	the	2-sigma	range and
% 95% range for	the	parameters S0, d and f in the ball and stick model in
% a single voxel

% Load and arrange the data
[dwis, qhat, bvals] = Load_DataSet_1();

% Select a voxel
Avox = dwis(:, 92, 65, 72);

% Define the number of bootstrap samples
T = 1000;

% Starting point
startx = [3.5e+00 3e-03 2.5e-01 0 0];

% Initialize the matrix where we are going to store the samples
total_samples = zeros(length(startx), T);

total_RESNORM = zeros(1, T);

% Run the bootstrap procedure
for t = 1 : T
    
    % Samples from the original measurements with replacement:
    indices = randi(length(Avox), length(Avox), 1);
    A_resampled = Avox(indices);
    qhat_resampled = qhat(:, indices);
    bvals_resampled = bvals(1, indices);

    % Inverse transformation to mantain the same starting point
    startx_n = startx;
    startx_n(1) = sqrt(startx(1));
    startx_n(2) = sqrt(startx(2));
    startx_n(3) = sqrt(-log(startx(3)));

    % Define various options for the non-linear fitting algorithm
    h=optimset('MaxFunEvals',20000,...
        'Algorithm','quasi-newton',...
        'MaxIter',2000,...
        'TolX',1e-10,...
        'TolFun',1e-10,...
        'Display','off');

    % Now run the fitting
    try
        [parameter_hat, total_RESNORM(1, t), ~, ~] = fminunc('BallStickSSD_constraints', startx_n, h, A_resampled, bvals_resampled, qhat_resampled);
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
sigma_ranges_cb = Two_Sigma_Range(total_samples);
ninety_five_range_cb = NinetyFive_Range(total_samples);
