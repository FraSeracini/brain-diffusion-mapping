function [parameters, dt, signal, Design_Matrix] = Linear_Diffusion_Tensor(Avox, bvals, qhat)

% This function calculates the diffusion tensor parameters based on a linear 
% estimation model, using log-transformed signal intensities and a least-squares 
% fitting approach.
%
% INPUTS:
%   - Avox : Vector of signal intensities for a given voxel
%   - bvals : Vector of b-values corresponding to each diffusion-weighted image
%   - qhat : Matrix of gradient directions (Nx3) associated with the diffusion images
%
% OUTPUTS:
%   - parameters : Estimated model parameters including log(S0) and tensor elements
%   - dt : The estimated 3x3 diffusion tensor matrix
%   - signal : Estimated signal from the fitted model
%   - Design_Matrix : The design matrix constructed from b-values and gradient directions
             

A = log(Avox);

% Construct the design matrix
Design_Matrix=[ones(1,length(bvals)); -bvals.*qhat(1,:).*qhat(1,:); -2*bvals.*qhat(1,:).*qhat(2,:); -2*bvals.*qhat(1,:).*qhat(3,:); -bvals.*qhat(2,:).*qhat(2,:); -2*bvals.*qhat(2,:).*qhat(3,:); -bvals.*qhat(3,:).*qhat(3,:) ]';

% Estimate parameters
parameters=Design_Matrix\A;

% Diffusion tensor
dt = [[parameters(2) parameters(3) parameters(4)]; [parameters(3) parameters(5) parameters(6)]; [parameters(4) parameters(6) parameters(7)]];

% Calculate the estimated signal
signal = exp(Design_Matrix * parameters);

end