function sumRes = BallStickSSD(x, Avox, bvals, qhat)
% This function implements the Ball-and-Stick diffusion model and calculates 
% the sum of squared differences between the measured diffusion-weighted signals 
% and the model-predicted signals for a single voxel.
%
% INPUTS:
%   - x : Vector of model parameters [S0, diff, f, theta, phi]
%       * S0 : Baseline signal without diffusion weighting
%       * diff : Diffusivity coefficient
%       * f : Volume fraction of the stick compartment
%       * theta : Orientation angle (elevation)
%       * phi : Orientation angle (azimuth)
%   - Avox : Measured signal intensities for the voxel
%   - bvals : Vector of b-values representing diffusion weightings
%   - qhat : Matrix of gradient directions (Nx3)
%
% OUTPUT:
%   - sumRes : Sum of squared differences (SSD) between measured and estimated signals

% Extract the parameters
S0 = x(1);
diff = x(2);
f = x(3);
theta = x(4);
phi = x(5);

% Synthesize the signals according to the model
fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat) 1])');
S = S0*(f*exp(-bvals*diff.*(fibdotgrad.^2)) + (1-f)*exp(-bvals*diff));

% Compute the sum of square differences
sumRes = sum((Avox - S').^2);
end