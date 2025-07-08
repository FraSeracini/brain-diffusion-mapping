%% V-Fold-Cross-Validation

% Load and arrange the new data set
[D, bvals, qhat] = Load_New_Dataset();

% Select a voxel
Avox = D(:, 1);

% Number of different starting points
starting_points = 20;

% Starting point
startx = [3.5e+00 3e-03 2.5e-01 0 0 3e-03];

% Select the number of V-folds
num_folds = 100;

% Calculate the V-Fold Cross-Validation errors
vfcv_errors = VFold_Cross_Validation('ZeppelinStickTortuosity', Avox, num_folds, startx, bvals, qhat, starting_points);

% Calculate the average error
average_vfcv_err = sum(vfcv_errors) / length(vfcv_errors)
