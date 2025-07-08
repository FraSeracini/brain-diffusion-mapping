%% Q1.3.3

% Use the AIC and the BIC to rank the four models

% Load and arrange the new data set
[D, bvals, qhat] = Load_New_Dataset();

% Select a voxel
Avox = D(:, 1);

% Load values
load('min_RESNORM_DT','SSD_DT');
load('min_RESNORM_BS','min_RESNORM_BS');
load('min_RESNORM_ZS','min_RESNORM_ZS');
load('min_RESNORM_ZST','min_RESNORM_ZST');

% Number of parameters in each model
n_parameters_DT = 7; 
n_parameters_BS = 5;
n_parameters_ZS = 6;
n_parameters_ZST = 5;

% Compute the AIC for each model
AIC_DT = Compute_AIC(n_parameters_DT, Avox, SSD_DT);
AIC_BS = Compute_AIC(n_parameters_BS, Avox, min_RESNORM_BS);
AIC_ZS = Compute_AIC(n_parameters_ZS, Avox, min_RESNORM_ZS);
AIC_ZST = Compute_AIC(n_parameters_ZST, Avox, min_RESNORM_ZST);

% Compute the AIC for each model
BIC_DT = Compute_BIC(n_parameters_DT, Avox, SSD_DT);
BIC_BS = Compute_BIC(n_parameters_BS, Avox, min_RESNORM_BS);
BIC_ZS = Compute_BIC(n_parameters_ZS, Avox, min_RESNORM_ZS);
BIC_ZST = Compute_BIC(n_parameters_ZST, Avox, min_RESNORM_ZST);