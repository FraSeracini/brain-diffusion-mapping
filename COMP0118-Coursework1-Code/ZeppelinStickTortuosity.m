function sumRes = ZeppelinStickTortuosity(x, Avox, bvals, qhat)
%
% This function fits the Zeppelin-and-Stick with Toruosity diffusion model to a single voxel's
% diffusion-weighted signals by minimizing the sum of squared differences between
% the measured and predicted signals.
%
% INPUTS:
%   - x : Vector of model parameters [S0, diff, f, theta, phi, s]
%       * S0 : Baseline signal intensity (constrained to be positive).
%       * diff : Diffusivity coefficient (constrained to be positive).
%       * f : Volume fraction of the stick compartment (constrained between [0,1]).
%       * theta : Elevation angle of fiber orientation.
%       * phi : Azimuth angle of fiber orientation.
%   - Avox : Measured signal intensities for the voxel.
%   - bvals : Vector of b-values representing diffusion weightings.
%   - qhat : Matrix of gradient directions (Nx3).

% Extract the parameters
S0 = (x(1)^2);
diff = (x(2)^2);
f = exp(-(x(3)^2));
theta = x(4);
phi = x(5);

% Synthesize the signals according to the model
fibdir = [cos(phi) * sin(theta) sin(phi) * sin(theta) cos(theta)];
fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat) 1])');
S = S0 * (f * exp(-bvals * diff .* (fibdotgrad.^2)) + (1 - f) .* exp(-bvals .* ((1-f).*diff+(diff-(1-f).*diff).* (fibdotgrad.^2))));

% Compute the sum
sumRes = sum((Avox - S').^2);

end