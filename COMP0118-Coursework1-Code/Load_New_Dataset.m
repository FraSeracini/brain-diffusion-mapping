function [D, bvals, qhat, TE] = Load_New_Dataset()
% LOAD_NEW_DATASET Switches to a different dataset and loads the diffusion signal and protocol.
%
% This function loads the diffusion signal data and the corresponding protocol parameters
% from 'isbi2015_data_normalised.txt' and 'isbi2015_protocol.txt'.
%
% OUTPUTS:
%   - meas : Column vector of diffusion signal measurements for the first voxel.
%   - bvals : Vector of computed b-values (in s/mm²).
%   - qhat : Matrix of gradient directions (3xN), where N is the number of measurements.

% Load the diffusion signal
fid = fopen('isbi2015_data_normalised.txt', 'r', 'b');
fgetl(fid); % Read in the header
D = fscanf(fid, '%f', [6, inf])'; % Read in the data
fclose(fid);

%% Load the protocol
fid = fopen('isbi2015_protocol.txt', 'r', 'b');
fgetl(fid);
A = fscanf(fid, '%f', [7, inf]);
fclose(fid);

% Create the protocol
grad_dirs = A(1:3, :);
G = A(4, :);
delta = A(5, :);
smalldel = A(6, :);
TE = A(7, :);
GAMMA = 2.675987E8;
bvals = ((GAMMA*smalldel.*G).^2).*(delta - smalldel / 3);

% Convert bvals units from s/m^2 to s/mm^2
bvals = bvals / 10^6;
qhat = A(1:3, :);

end