function FIM = Fisher_Information_Matrix(parameters, Avox, bvals, qhat)
%
% This function calculates the Fisher Information Matrix (FIM) based on the 
% sensitivity of the predicted signal to model parameters, providing a measure 
% of parameter uncertainty.
%
% INPUTS:
%   - parameters : Estimated model parameters [S0, d, f, theta, phi].
%   - Avox       : Measured diffusion-weighted signal intensities.
%   - bvals      : Vector of b-values representing diffusion weightings.
%   - qhat       : Matrix of gradient directions (Nx3).
%
% OUTPUT:
%   - FIM : Fisher Information Matrix (num_params x num_params), where 
%           num_params corresponds to the estimated model parameters.

% Considering only the first 3 parameters (ignore orientation parameters)
num_params = 3;

% Initialize the FIM
FIM = zeros(num_params, num_params);

% % Conversione dell'unità di diffusione in 10⁻³ mm²/s per migliorare la stabilità numerica
parameters(2) = parameters(2) * 1e3;

% Synthetize the estimated signal
A_est = BallStickSSD_constraints_signal(parameters, bvals, qhat);

% Estimate the variance
sigma2 = (1 / (length(Avox) - num_params)) * sum((Avox - A_est').^2);


% Delta for the numerical calculation of the derivative
epsilon = 1e-3 * abs(parameters);

% Initialize the matrix where we are going to store the results of the
% partial derivatives
dSdx = zeros(length(A_est), num_params);

% Calculate all the numerical derivatives
for i = 1:num_params
    parameters_increased = parameters;
    parameters_increased(i) = parameters_increased(i) + epsilon(i);
    A_increased = BallStickSSD_constraints_signal(parameters_increased, bvals, qhat);
    dSdx(:, i) = (A_increased - A_est) / epsilon(i);
end


% Normalize dSdx to avoid numerical instability
dSdx = dSdx ./ std(dSdx, [], 1);
% Compute the FIM
for i = 1:num_params
    for j = i:num_params
        FIM(i, j) = sum((1 / sigma2) * dSdx(:, i) .* dSdx(:, j));
        FIM(j, i) = FIM(i, j); % Assure the symmetry of the matrix
    end
end

% Normalize FIM using xxᵀ (J^T * J)
FIM = (dSdx' * dSdx) / sigma2;


end
