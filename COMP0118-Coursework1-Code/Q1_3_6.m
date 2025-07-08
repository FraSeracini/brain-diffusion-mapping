%% Q1.3.6
% Fitting the simple two compartment Zeppelin-and-stick	models in Q1.3.2 to
% the imaging data we used in sections 1.1 and 1.2. Use the AIC to select
% the most appropriate model in each voxel of the image	and	create a map
% indicating which model is most appropriate at each location

% Load and arrange the data
[dwis, qhat, bvals] = Load_DataSet_1();

% Number of different starting points for each optimization process
starting_points = 15;

% Initialize the colour maps
colour_map_AIC = zeros(145, 174, 3);
colour_map_BIC = zeros(145, 174, 3);
% colours = [0 1 0; 0 0 1; 1 0 0];
colours = [0 0 1; 1 0 0];

% Initialize the matrix
best_parameters_BS = zeros(145, 174, 5);
best_parameters_ZS = zeros(145, 174, 6);
best_parameters_ZST = zeros(145, 174, 5);
min_RESNORM_ZST = zeros(145, 174);
min_RESNORM_ZS = zeros(145, 174);
min_RESNORM_BS = zeros(145, 174);
AIC_BS = zeros (145, 174);
% AIC_ZS = zeros (145, 174);
% AIC_ZST = zeros (145, 174);
% BIC_BS = zeros (145, 174);
% BIC_ZS = zeros (145, 174);
% BIC_ZST = zeros (145, 174);


% Map the AIC for the different models in each voxel of the chosen slice
for i = 1 : 145
    for j = 1 : 174

        % Select the voxel
        Avox = dwis(:, i, j, 72);

        % Use the results from the dt to create a starting point
        % startx_ZST = DT_starting_point(Avox, bvals, qhat);
        % startx_ZS = [startx_ZST startx_ZST(3)];
        % startx_BS = startx_ZST;

        if min(Avox) > 0 % Accept only significant values

            % Run the fitting
            % [best_parameters_ZST(i, j, :), min_RESNORM_ZST(i, j), ~] = find_optimal_parameters('ZeppelinStickTortuosity', startx_ZST, starting_points, Avox, bvals, qhat);
            % [best_parameters_ZS(i, j, :), min_RESNORM_ZS(i, j), ~] = optimal_fit_ZeppelinStick_first_dataset(startx_ZS, starting_points, Avox, bvals, qhat);
            % [best_parameters_BS(i, j, :), min_RESNORM_BS(i, j), ~] = find_optimal_parameters('BallStickSSD_constraints', startx_BS, starting_points, Avox, bvals, qhat);

            % % Compute the AIC for the different models
            % AIC_BS(i, j) = Compute_AIC(5, Avox, min_RESNORM_BS(i, j));
            % AIC_ZS(i, j) = Compute_AIC(6, Avox, min_RESNORM_ZS(i, j));
            % AIC_ZST(i, j) = Compute_AIC(5, Avox, min_RESNORM_ZST(i, j));
            % 
            % % Compute the BIC for the different models
            % BIC_BS(i, j) = Compute_BIC(5, Avox, min_RESNORM_BS(i, j));
            % BIC_ZS(i, j) = Compute_BIC(6, Avox, min_RESNORM_ZS(i, j));
            % BIC_ZST(i, j) = Compute_BIC(5, Avox, min_RESNORM_ZST(i, j));

            % Assign colours to the AIC colour map
            % [~, index] = min([AIC_BS(i, j), AIC_ZS(i, j), AIC_ZST(i, j)]);
            % [~, index] = min([AIC_ZS(i, j), AIC_ZST(i, j)]);
            % colour_map_AIC(i, j, :) = colours(index, :);

            if BIC_ZS(i,j) < BIC_ZST(i,j)
                colour_map_BIC(i,j,:) = [0 0 1];
            elseif BIC_ZST(i,j) < BIC_ZS(i,j)
                colour_map_BIC(i,j, :) = [1 0 0];
            end
            
            % Assign colours to the BIC colour map
            % [~, index] = min([BIC_BS(i, j), BIC_ZS(i, j), BIC_ZST(i, j)]);
            % [~, index] = min([BIC_ZS(i, j), BIC_ZST(i, j)]);
            % colour_map_BIC(i, j, :) = colours(index, :);
            

        end
    end

end

figure
imshow(flipud(permute(colour_map_AIC, [2, 1, 3])), [])
title('AIC Map')

figure
imshow(flipud(permute(colour_map_BIC, [2, 1, 3])), [])
title('BIC Map')

%% Average parameters estimates from different models using Akaike weights

% Initialize the matrix where we are going to store the weighted RESNORM
% values
RESNORM_map = zeros(145, 174);

% Initialize the matrix where we are going to store the weighted mean of
% the parameters
best_param_mean = zeros(145, 174, 5);

% Apply the Akaike weights to all the voxels in the chosen slice
for i = 1 : 145
    for j = 1 : 174
        
        % Select a voxel
        Avox = dwis(:, i ,j, 72);

        % Calculate the minimun AIC between the 3 models
        AIC_min = min([AIC_BS(i, j), AIC_ZS(i, j), AIC_ZST(i, j)]);
       
        % Calculate the AIC differences
        AIC_diff_BS = AIC_BS(i, j) - AIC_min;
        AIC_diff_ZS = AIC_BS(i, j) - AIC_min;
        AIC_diff_ZST = AIC_BS(i, j) - AIC_min;

        % Calculate the Akaike weights
        Ak_BS = exp(-0.5 * AIC_diff_BS) / (exp(-0.5 * AIC_diff_BS) + exp(-0.5 * AIC_diff_ZS)+ exp(-0.5 * AIC_diff_ZST));
        Ak_ZS = exp(-0.5 * AIC_diff_ZS) / (exp(-0.5 * AIC_diff_BS) + exp(-0.5 * AIC_diff_ZS)+ exp(-0.5 * AIC_diff_ZST));
        Ak_ZST = exp(-0.5 * AIC_diff_ZST) / (exp(-0.5 * AIC_diff_BS) + exp(-0.5 * AIC_diff_ZS)+ exp(-0.5 * AIC_diff_ZST));

        % Calculate the weighted mean of the optimal parameters
        best_param_mean(i, j, :) = Ak_BS * best_parameters_BS(i, j, :) + Ak_ZS * best_parameters_ZS(i, j, 1:5) + Ak_ZST * best_parameters_ZST(i, j, :);

        
        % Calculate the estimated signal for the 3 models using the
        % weighted mean of the optimal parameters
        Aest_BS = BallStickSSD_constraints_signal(best_param_mean(i, j, :), bvals, qhat);
        Aest_ZS = ZeppelinStickSSD_signal([squeeze(best_param_mean(i, j, :))', best_parameters_ZS(i, j, 6)], bvals, qhat);
        Aest_ZST = ZeppelinStickTortuosity_signal(best_param_mean(i, j, :), bvals, qhat);

        % Calculate the weighted mean of the estimated signals
        Aest_weighted = Ak_BS .* Aest_BS + Ak_ZS .* Aest_ZS + Ak_ZST .* Aest_ZST;

        % Calculate the sum and store it in the RESNORM map
        RESNORM_map(i, j) = sum((Avox - Aest_weighted').^2);

    end
end

% Visualize the RESNORM map
figure
imshow(flipud(RESNORM_map'), [])
title('RESNORM map')