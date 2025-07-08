function ninety_five_range = NinetyFive_Range(total_samples)
%
% This function calculates the 95% confidence interval for the parameters 
% S0, d, and f by sorting the sample values and selecting the 2.5th and 
% 97.5th percentiles.
%
% INPUT:
%   - total_samples : A matrix where each row contains sampled values of 
%                     a specific parameter (S0, d, or f) from multiple runs.
%
% OUTPUT:
%   - ninety_five_range : A 3x2 matrix containing the lower and upper bounds 
%                         of the 95% confidence interval for each parameter.

% Initialize the matrix where we are going to store the 95% range
ninety_five_range = zeros(3, 2);

% Sort in ascending order the samples
sorted_samples = sort(total_samples, 2);

% Estimate the 95% range for the parameters S0, d and f
for i = 1 : 3
    sorted_samples_no_zeros = sorted_samples(i, sorted_samples(i, :) ~= 0);
    ninety_five_range(i, 1) = sorted_samples(i, round(length(sorted_samples_no_zeros) * 0.025));
    ninety_five_range(i, 2) = sorted_samples(i, round(length(sorted_samples_no_zeros) * (1 - 0.025)));
end
end