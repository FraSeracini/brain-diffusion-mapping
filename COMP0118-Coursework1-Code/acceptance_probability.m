function alpha = acceptance_probability(sigma, Avox, bvals, qhat, Xc, Xt)
% This function calculates the acceptance probability used in the Metropolis-Hastings 
% step of the Markov Chain Monte Carlo (MCMC) algorithm for fitting the Ball-and-Stick model.
%
% INPUTS:
%   - sigma : Standard deviation of the noise model.
%   - Avox : Measured signal intensities for the voxel.
%   - bvals : Vector of b-values representing diffusion weightings.
%   - qhat : Matrix of gradient directions (Nx3).
%   - Xc : Proposed new parameter set [S0, d, f, theta, phi].
%   - Xt : Current parameter set [S0, d, f, theta, phi].
%
% OUTPUT:
%   - alpha : Acceptance probability for the MCMC update step.

% Extract the estimated signals for both set of parameters
Signal_1 = BallStickSSD_constraints_signal(Xc, bvals, qhat);
Signal_2 = BallStickSSD_constraints_signal(Xt, bvals, qhat);

%Extract the theta parameters 
theta_c = Xc(4);
theta_t = Xt(4);

% Calculate the acceptance probabilty
alpha = min(1, (sin(theta_c) / sin(theta_t)) * exp(-(sum((Avox - Signal_1').^2) - sum((Avox - Signal_2').^2)) / (2 * sigma^2)));

end