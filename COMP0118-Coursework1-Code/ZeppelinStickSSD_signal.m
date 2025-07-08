function S = ZeppelinStickSSD_signal(x, bvals, qhat)
%
% This function synthesizes diffusion-weighted signals based on the 
% Zeppelin-and-Stick model parameters for a single voxel.
%
% INPUTS:
%   - x : Vector of model parameters [S0, diff, f, theta, phi, s]
%       * S0 : Baseline signal intensity.
%       * diff : Diffusivity coefficient.
%       * f : Volume fraction of the stick compartment.
%       * theta : Elevation angle of fiber orientation.
%       * phi : Azimuth angle of fiber orientation.
%       * s : Scaling factor to calculate the smallest eigenvalue.
%   - bvals : Vector of b-values representing diffusion weightings.
%   - qhat : Matrix of gradient directions (Nx3).
%
% OUTPUT:
%   - S : Vector of predicted signal intensities for each diffusion-weighted measurement.

% Extract the parameters
S0 = x(1);
diff = (x(2));
f = x(3);
theta = x(4);
phi = x(5);
s = x(6); 
lambda2 = s * diff; 


% Synthesize the signals according to the model
fibdir = [cos(phi) * sin(theta) sin(phi) * sin(theta) cos(theta)];
fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat) 1])');
S = S0 * (f * exp(-bvals .* diff .* (fibdotgrad.^2)) + (1 - f) * exp(-bvals .* (lambda2 + (diff - lambda2) .* (fibdotgrad.^2))));

end