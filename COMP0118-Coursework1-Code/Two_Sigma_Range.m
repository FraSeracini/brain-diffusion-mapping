function sigma2_ranges = Two_Sigma_Range(total_samples)
%
% This function calculates the 2-sigma range (95.4% confidence interval) for 
% the parameters S0, d, and f, by removing zero values, computing the mean, 
% and applying the standard deviation.
%
% INPUT:
%   - total_samples : A matrix where each row contains sampled values of 
%                     a specific parameter (S0, d, or f) from multiple runs.
%
% OUTPUT:
%   - sigma2_ranges : A 3x2 matrix containing the lower and upper bounds 
%                     of the 2-sigma range for each parameter.

% Initialize the matrix where we are going to store the 2-sigma range
sigma2_ranges = zeros(3, 2);

% Estimate the 2-sigma range for the parameters S0, d and f
for i = 1 : 3
    total_samples_no_zeros = total_samples(i, total_samples(i, :) ~= 0);
    sigma2_ranges(i, 1) = mean(total_samples_no_zeros) - 2 * std(total_samples_no_zeros);
    sigma2_ranges(i, 2) = mean(total_samples_no_zeros) + 2 * std(total_samples_no_zeros);
end
end