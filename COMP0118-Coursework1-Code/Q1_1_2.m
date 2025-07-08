%% Q1.1.2

% Ball-and-Stick Model to measure, indirectly, the density of axon fibres in
% each voxel and thus produce a map	of that parameter over the brain.

% Load and arrange the data
[dwis, qhat, bvals] = Load_DataSet_1();

% Select a voxel
Avox = dwis(:, 92, 65, 72);

% Define a starting point for the non-linear fit
startx = [3.5e+00   3e-03   2.5e-01 0 0];

% Define various options for the non-linear fitting algorithm.
h=optimset('MaxFunEvals',10000,...
    'Algorithm','quasi-newton',...
    'MaxIter', 2000,...
    'TolX',1e-10,...
    'TolFun',1e-10,...
    'Display','Iter');

% Now run the fitting
[parameter_hat, RESNORM, ~, ~] = fminunc('BallStickSSD', startx, h, Avox, bvals, qhat);

% Visualize the results
parameter_hat
RESNORM

% Extract the estimated signal with the parameters obtained
Aest = BallStickSSD_constraints_signal(parameter_hat, bvals, qhat);

% Visualize the fit to the data by comparing the predictions from the model
% with the fitted parameters obtained to the actual measurements
figure;
plot(Avox, 'bs')
hold on;
plot(Aest, 'rx')
xlabel('k');
ylabel('S');
legend('Data', 'Model');
