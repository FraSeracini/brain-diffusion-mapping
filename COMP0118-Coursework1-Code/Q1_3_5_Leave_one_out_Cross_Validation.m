%% Leave-one-out-Cross-Validation

% Load and arrange the new data set
[D, bvals, qhat] = Load_New_Dataset();

% Select a voxel
Avox = D(:, 1);

% Number of different starting points
starting_points = 20;

% Starting point
startx = [3.5e+00 3e-03 2.5e-01 0 0 3e-03];

% Calculate the Leave-one-out Cross-Validation errors
l1ocv_errors = L1_out_Cross_Validation('ZeppelinStickSSD', Avox, startx, bvals, qhat, starting_points);

% Calculate the average error
average_l1ocv = sum(l1ocv_errors) / length(l1ocv_errors)