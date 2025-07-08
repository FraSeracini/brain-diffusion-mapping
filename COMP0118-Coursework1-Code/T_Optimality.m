function T_optimality = T_Optimality(parameters, Avox, bvals, qhat)
%
% This function calculates the T-optimality value, which is defined as the 
% trace of the Fisher Information Matrix (FIM). A higher value 
% of T-optimality indicates more precise parameter estimates.
%
% INPUTS:
%   - parameters : Estimated model parameters.
%   - Avox       : Experimental diffusion-weighted signal data.
%   - bvals      : Vector of b-values representing diffusion weightings.
%   - qhat       : Matrix of gradient directions (Nx3).
%
% OUTPUT:
%   - T_optimality : T-optimality value (trace of FIM).

% Calculate the FIM
FIM = Fisher_Information_Matrix(parameters, Avox, bvals, qhat);

% Calculate the T-optimality value
T_optimality = trace(FIM);

end
