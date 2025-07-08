%% Q1.1.3

% Adapt the implementation in Q1.1.2 to use the transformation method
% to allow only physically realistic settings of the parameters: S0 must be
% positive, d must be positive, and f must be in [0, 1]

% Load and arrange the data
[dwis, qhat, bvals] = Load_DataSet_1();

% Select a voxel
Avox = dwis(:, 92, 65, 72);

% Define a starting point for the non-linear fit
startx_constraints = [sqrt(3.5e+00) sqrt(3e-03) sqrt(-log(0.001e-01)) 0 0];

% Define various options for the non-linear fitting algorithm.
h=optimset('MaxFunEvals',10000,...
    'Algorithm','quasi-newton',...
    'MaxIter', 2000,...
    'TolX',1e-10,...
    'TolFun',1e-10,...
    'Display','off');

% Now run the fitting
[parameter_hat, RESNORM, EXITFLAG, OUTPUT] = fminunc('BallStickSSD_constraints', startx_constraints, h, Avox, bvals, qhat);

% Transformation to get the model parameters
parameter_hat(1) = parameter_hat(1)^2;
parameter_hat(2) = parameter_hat(2)^2;
parameter_hat(3) = exp(-parameter_hat(3)^2);

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
