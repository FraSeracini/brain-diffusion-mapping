function S = ZeppelinStickTortuosity_signal(x, bvals, qhat)
%
% This function synthesizes diffusion-weighted signals based on the 
% Zeppelin-and-Stick with Tortuosity model parameters for a single voxel.
%
% INPUTS:
%   - x : Vector of model parameters [S0, diff, f, theta, phi]
%       * S0 : Baseline signal intensity.
%       * diff : Diffusivity coefficient.
%       * f : Volume fraction of the stick compartment.
%       * theta : Elevation angle of fiber orientation.
%       * phi : Azimuth angle of fiber orientation.
%   - bvals : Vector of b-values representing diffusion weightings.
%   - qhat : Matrix of gradient directions (Nx3).
%
% OUTPUT:
%   - S : Vector of predicted signal intensities for each diffusion-weighted measurement.

% Extract the parameters
S0 = x(1);
diff = x(2);
f = x(3);
theta = x(4);
phi = x(5);

% Synthesize the signals according to the model
fibdir = [cos(phi) * sin(theta) sin(phi) * sin(theta) cos(theta)];
fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat) 1])');
S = S0 * (f * exp(-bvals * diff .* (fibdotgrad.^2)) + (1 - f) .* exp(-bvals .* ((1-f).*diff+(diff-(1-f).*diff).* (fibdotgrad.^2))));

end