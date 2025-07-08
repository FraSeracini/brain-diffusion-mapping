function BIC = Compute_BIC(n_parameters, Avox, SSD)
%
% This function calculates the BIC value to assess the goodness of fit of a model 
% while penalizing model complexity.
%
% INPUTS:
%   - n_parameters : Number of estimated model parameters.
%   - Avox : Measured signal intensities for the voxel.
%   - SSD : Sum of squared differences (SSD) between measured and predicted signals.
%
% OUTPUT:
%   - BIC : Computed Bayesan Information Criterion value

% Accounting for the estimate of sigma
n_parameters = n_parameters+1;

% Calculate the value of K for the AIC
K = length(Avox);

% Use an approsimation to calculate sigma^2
sigma2 = (1 / K) * SSD;

% Compute the BIC
BIC = n_parameters * log(K) + K * log(sigma2);

end