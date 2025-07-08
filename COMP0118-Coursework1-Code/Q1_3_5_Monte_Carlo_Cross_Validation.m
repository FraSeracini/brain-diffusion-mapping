%% Monte-Carlo Cross-Validation

% Load and arrange the new data set
[D, bvals, qhat] = Load_New_Dataset();

% Select a voxel
Avox = D(:, 1);

% Number of different starting points
starting_points = 20;

% Starting point
startx = [3.5e+00 3e-03 2.5e-01 0 0 3e-03];

% Select the number of iterations
num_iter = 500;

% Select the dimension of the test set
test_dimension = 20;

% Calculate the Monte-Carlo Cross-Validation errors
mccv_errors = MCarlo_Cross_Validation('ZeppelinStickSSD', Avox, test_dimension, num_iter, startx, bvals, qhat, starting_points);

%Calculate the average error
average_mccv_err = sum(mccv_errors) / length(mccv_errors)