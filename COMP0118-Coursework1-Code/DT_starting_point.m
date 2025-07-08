function startx = DT_starting_point(Avox, bvals, qhat)
%
% This function estimates the diffusion tensor for a given voxel and extracts 
% key parameters to serve as a starting point for non-linear optimization in 
% models such as the Ball-and-Stick model.
%
% INPUTS:
%   - Avox : Measured signal intensities for the voxel.
%   - bvals : Vector of b-values representing diffusion weightings.
%   - qhat : Matrix of gradient directions (Nx3).
%
% OUTPUT:
%   - startx : Initial parameter estimates for the Ball-and-Stick model.
%              [S0, mean diffusivity, volume fraction, theta, phi]

% Estimate the Linear Diffusion Tensor model
[x, dt, ~, ~] = Linear_Diffusion_Tensor(Avox, bvals, qhat);

% Find the eigenvalues of the Diffusion Tensor Matrix
[eigVec, eigVal] = eig(dt);

% Find the principal direction of the diffusivity
[~, index] = max(diag(eigVal));
direction = eigVec(:, index);

% Transform the main direction in spherical angles
theta = acos(direction(3));
phi = atan2(direction(2), direction(1));

% Use the results from the dt to create a starting point for the
% ball-and-stick model
startx = [exp(x(1)) trace(dt)/3 2.5e-01 theta phi];

end