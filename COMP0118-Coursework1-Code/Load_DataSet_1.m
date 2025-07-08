function [dwis, qhat, bvals] = Load_DataSet_1()
% This function loads the diffusion-weighted imaging (DWI) data, gradient 
% directions, and b-values for further processing in diffusion modeling.
%
% OUTPUTS:
%   - dwis : 4D matrix containing diffusion-weighted imaging data, 
%            rearranged for proper indexing.
%   - qhat : Matrix of gradient directions used in the MRI acquisition.
%   - bvals : Vector of b-values computed from gradient directions, 
%             indicating the diffusion weighting applied in each measurement.
%
% The function assumes the presence of 'data' and 'bvecs' files in the 
% working directory.

% LOAD AND ARRANGE THE DATA
load('data');
dwis = double(dwis);
dwis = permute(dwis,[4, 1, 2, 3]);
qhat = load('bvecs');
bvals = 1000*sum(qhat.*qhat);

end
