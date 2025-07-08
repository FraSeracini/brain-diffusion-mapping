function sumRes = BallStickSSD_constraints(x, Avox, bvals, qhat)
% This function implements a constrained version of the Ball-and-Stick diffusion model,
% ensuring physically meaningful parameter values. It calculates the sum of squared
% differences between the measured and model-predicted diffusion-weighted signals
% for a single voxel.
%
% INPUTS:
%   - x : Vector of constrained model parameters [S0, diff, f, theta, phi]
%       * S0 : Baseline signal intensity (constrained to be positive)
%       * diff : Diffusivity coefficient (constrained to be positive)
%       * f : Volume fraction of the stick compartment (constrained between [0,1])
%       * theta : Orientation angle (elevation)
%       * phi : Orientation angle (azimuth)
%   - Avox : Measured signal intensities for the voxel
%   - bvals : Vector of b-values representing diffusion weightings
%   - qhat : Matrix of gradient directions (Nx3)
%
% OUTPUT:
%   - sumRes : Sum of squared differences (SSD) between measured and predicted signal

% Extract the parameters
S0 = (x(1)^2);
diff = (x(2)^2);
f = exp(-(x(3)^2));
theta = x(4);
phi = x(5);

% Synthesize the signals according to the model
fibdir = [cos(phi) * sin(theta) sin(phi) * sin(theta) cos(theta)];
fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat) 1])');
S = S0 * (f * exp(-bvals * diff .* (fibdotgrad.^2)) + (1 - f) * exp(-bvals * diff));

% Compute the sum
sumRes = sum((Avox - S').^2);
end