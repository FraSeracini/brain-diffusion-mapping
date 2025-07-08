function A_optimality = A_Optimality(parameters, Avox, bvals, qhat)
%
% This function calculates the A-optimality value, which is defined as the 
% trace of the inverse of the Fisher Information Matrix (FIM). A lower value 
% of A-optimality indicates more precise parameter estimates.
%
% INPUTS:
%   - parameters : Estimated model parameters.
%   - Avox       : Experimental diffusion-weighted signal data.
%   - bvals      : Vector of b-values representing diffusion weightings.
%   - qhat       : Matrix of gradient directions (Nx3).
%
% OUTPUT:
%   - A_optimality : A-optimality value (trace of FIM^-1).

% Calculate the FIM
FIM = Fisher_Information_Matrix(parameters, Avox, bvals, qhat);

% Calculate the A-optimality value
A_optimality = trace(inv(FIM));

end
