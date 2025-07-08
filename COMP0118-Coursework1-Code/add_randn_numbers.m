function x_n = add_randn_numbers(x)
% For the data set used in sections 1.1 and 1.2
%
% This function perturbs the input parameter vector by adding normally 
% distributed random noise to generate new initial conditions for optimization.
% It is used to explore different starting points in the optimization process 
% for diffusion model fitting.
%
% INPUT:
%   - x : Original parameter vector [S0, d, f, theta, phi]
%
% OUTPUT:
%   - x_n : New perturbed parameter vector with added random noise.

% Perturb the parameters
x_n(1) = x(1) + randn*1e+03;
x_n(2) = x(2) + randn*1e-03;
x_n(3) = x(3) + randn*1e-01;
x_n(4) = x(4) + randn*1e-01;
x_n(5) = x(5) + randn*1e-01;

end