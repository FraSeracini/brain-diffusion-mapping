function shells = Find_Shells(bvals, TE)
%
% This function categorizes diffusion MRI acquisitions into unique shells 
% based on their b-values and groups them with corresponding b=0 acquisitions 
% that share the same echo time (TE).
%
% INPUTS:
%   - bvals : Vector of b-values for diffusion acquisitions.
%   - TE    : Vector of echo times associated with each acquisition.
%
% OUTPUT:
%   - shells : A structured array where each entry corresponds to a shell 
%              and contains:
%              * b_value : The unique b-value of the shell.
%              * indices : Indices of acquisitions belonging to this shell, 
%                          including corresponding b=0 scans.

% Find unique b-values (temporarily excluding b=0)
all_bvals = unique(bvals(bvals > 0));

% Initialize a structure for the shells
shells = struct();

% Find acquisitions with b = 0
b0_index = find(bvals == 0);

for i = 1:length(all_bvals)
    % Find measurements with this b-value
    shell_index = find(bvals == all_bvals(i));

    % Find TE values for this shell
    TE_shell = unique(TE(shell_index));

    % Select b = 0 acquisitions with the same TE
    b0_matching_index = b0_index(ismember(TE(b0_index), TE_shell));

    % Save the shell indices
    shells(i).b_value = all_bvals(i);
    shells(i).indices = [shell_index b0_matching_index];
end
end
