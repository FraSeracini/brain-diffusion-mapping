%% Q1.1.1

% Implement Linear Diffusion Tensor Estimator and use it to map the mean
% diffusivity, fractional anisotropy, and colour-coded principal direction
% over one slice of the image.

% Load and arrange the data
[dwis, qhat, bvals] = Load_DataSet_1();

% Put data for one voxel into A
Avox = dwis(:, 92, 65, 72);

% Calculate the parameters for the diffusion tensor
[x, ~, ~, Y] = Linear_Diffusion_Tensor(Avox, bvals, qhat);

% For all the voxels in slice 72

% Initialize the maps
dt_map72 = zeros(7, 145, 174);
fa = zeros(145, 174);
colour_map = zeros(145, 174, 3);

% Compute the Mean Diffusivity Map
for i = 1 : 145
    for j = 1 : 174
        A = dwis(:, i, j, 72);

        if(min(A) >0) % Accept only significant values

            % Diagonal matrix of weights
            W = diag(A.^2);

            % Compute the inverse map
            invmap = inv(Y' * W * Y) * Y' * W;

            % Estimate parameters
            x = invmap * log(A);

            % Save parameters in the Mean Diffusivity Map
            dt_map72(:, i, j) = x;

        end
    end
end

% Compute the Fractional Anisotropy Map and the colour-coded principal
% direction map
for i = 1 : 145
    for j = 1 : 174

        % Retrieve the parameters from the Mean Diffusivity Map
        x = dt_map72(:, i, j);

        % Save the parameters in the Diffusion Tensor
        dt = [[x(2) x(3) x(4)]; [x(3) x(5) x(6)]; [x(4) x(6) x(7)]];

        % Find the eigenvalues of the Diffusion Tensor Matrix
        [eigVec, eigVal] = eig(dt);

        % Compute the Fractional Anisotropy Map
        fa(i, j) = sqrt(1.5 * sum((diag(eigVal) - trace(dt) / 3).^2) / sum(diag(eigVal).^2));

        % Find the principal direction of the diffusivity
        [~, index] = max(diag(eigVal));
        direction = eigVec(:, index);

        % Compute the colour-coded principal direction map
        colour_map(i, j, :) = (fa(i, j)^0.5) * abs(direction);
    end
end

% Visualize the maps
figure;
imshow(flipud(squeeze(dt_map72(2, :, :)  + dt_map72(5, :, :) + dt_map72(7, :, :))'),[])
title('Mean Diffusivity Map')

figure;
imshow(flipud(fa'))
title('Fractional Anisotropy Map')

figure;
imshow(flipud(pagetranspose(colour_map)))
title('Directionally encoded colour map')
toc
