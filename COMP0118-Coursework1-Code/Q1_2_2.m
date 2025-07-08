%% Q1.2.2

% Use MCMC to provide another estimate of the 2-sigma range and the
% 95% range for the same parameters

% Load and arrange the data
[dwis, qhat, bvals] = Load_DataSet_1();

% Select a voxel
Avox = dwis(:, 92, 65, 72);

% Starting points
starting_points = 15;

% Select the number of samples
T = 1000000;

% Select the number of burn-in samples
burn_in = 20000;

% Select the sampling interval
Sampling_interval = 2000;

% Starting point
startx = [3.5e+00 3e-03 2.5e-01 0 0];

% Run the fitting procedure from different starting points
if(min(Avox) > 0) % Accept only significant values

    % Find the parameters associated to global minima of the
    % objective function
    [startx, ~, ~] = find_optimal_parameters('BallStickSSD_constraints', startx, starting_points, Avox, bvals, qhat);

end

% Extract the estimated signal with the parameters obtained
A_est = BallStickSSD_constraints_signal(startx, bvals, qhat);

% Compute the standard deviation
sigma = sqrt((1 / (length(Avox) - length(startx))) * sum((Avox - A_est').^2));

% Initialize the matrix where we are going to store the samples
mcmc_samples = zeros(T, length(startx));
mcmc_samples(1, :) = startx;

% Apply the MCMC algorithm
for t = 2 : T

    % Random points drawn from proposal distribution
    Xc(1) = normrnd(mcmc_samples(t-1, 1), 1e+02);
    Xc(2) = normrnd(mcmc_samples(t-1, 2), 1e-05);
    Xc(3) = normrnd(mcmc_samples(t-1, 3), 1e-05);
    Xc(4) = normrnd(mcmc_samples(t-1, 4), 1e-02);
    Xc(5) = normrnd(mcmc_samples(t-1, 5), 1e-02);

    % Impose that S0 >= 0
    while (Xc(1) < 0)
        Xc(1) = normrnd(mcmc_samples(t-1, 1),1e+01);
    end

    % Impose that d >= 0
    while (Xc(2) < 0)
        Xc(2) = normrnd(mcmc_samples(t-1,2),1e-05);
    end

    % Impose that 0 <= f <= 1
    while (Xc(3) < 0 || Xc(3) >1)
        Xc(3) = normrnd(mcmc_samples(t-1,3),1e-05);
    end
    
    % Calculate the acceptance probability
    alpha = acceptance_probability(sigma, Avox, bvals, qhat, Xc, mcmc_samples(t-1, :));
    
    % Accept the sample only if the acceptance probability is greater than a
    % random probability drawn from a uniform distribution
    if alpha > rand
        mcmc_samples(t, :) = Xc;
    else
        mcmc_samples(t, :) = mcmc_samples(t-1, :);
    end
end

% Initialize the matrix where we are going to store only the interesting
% samples
final_mcmc_samples = zeros((T-burn_in) / Sampling_interval, length(startx));

i = 1;
% Discard the first burn-in samples and accept only the samples multiples
% of the sampling interval
for t = burn_in : Sampling_interval : (T - burn_in)
    final_mcmc_samples(i, :) = mcmc_samples(t, :);
    i = i + 1;
end

% Transpose the samples
final_mcmc_samples = final_mcmc_samples';

% Estimate the 2-sigma range and the 95% range for the parameters S0, d and f
sigma_ranges_mcmc = Two_Sigma_Range(final_mcmc_samples);
ninety_five_range_mcmc = NinetyFive_Range(final_mcmc_samples);
