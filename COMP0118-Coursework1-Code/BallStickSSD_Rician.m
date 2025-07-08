function sumRes = BallStickSSD_Rician(x, Avox, bvals, qhat)
%
% This function fits the Ball-and-Stick diffusion model while accounting for 
% the presence of Rician noise using the Corrected Least Squares (CLS) approach.
%
% INPUTS:
%   - x : Vector of model parameters [S0, diff, f, theta, phi, sigma]
%       * S0 : Baseline signal intensity.
%       * diff : Diffusivity coefficient.
%       * f : Volume fraction of the stick compartment.
%       * theta : Elevation angle of fiber orientation.
%       * phi : Azimuth angle of fiber orientation.
%       * sigma : Standard deviation of Rician noise (ensured to be positive).
%   - Avox : Measured diffusion-weighted signal intensities.
%   - bvals : Vector of b-values representing diffusion weightings.
%   - qhat : Matrix of gradient directions (Nx3).
%
% OUTPUT:
%   - sumRes : Sum of squared differences (SSD) between the measured signals 
%              and the predicted signals corrected for Rician noise.

% Extract the parameters
S0 = x(1);
diff = x(2);
f = x(3);
theta = x(4);
phi = x(5);
sigma = abs(x(6)); % Ensure positivity of noise parameter

% Synthesize the signals according to the model
fibdir = [cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)];
fibdotgrad = sum(qhat .* repmat(fibdir, [length(qhat), 1])', 1);
S = S0 * (f * exp(-bvals * diff .* (fibdotgrad.^2)) + (1 - f) * exp(-bvals * diff));

% Apply Rician noise correction
S_corrected = sqrt(S.^2 + sigma^2);

% Compute the sum 
sumRes = sum((Avox - S_corrected').^2);

end
