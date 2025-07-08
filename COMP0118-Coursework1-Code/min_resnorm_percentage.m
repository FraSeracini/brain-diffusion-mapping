function [percentage, num_runs] = min_resnorm_percentage(RESNORM_per_start_point, starting_points, accepted_deviation)
%
% This function calculates the proportion of trials where the smallest RESNORM 
% value is found and estimates the required number of runs to achieve a 95% 
% probability of finding the global minimum at least once.
%
% INPUTS:
%   - RESNORM_per_start_point : Array containing RESNORM values from each optimization run.
%   - starting_points : Total number of optimization runs performed.
%   - accepted_deviation : Maximum allowed deviation from the global minimum.
%
% OUTPUTS:
%   - percentage : Proportion of trials that found the minimum RESNORM value.
%   - num_runs : Estimated number of runs required to have a 95% probability of 
%                finding the global minimum at least once.

% Find the global minima
min_RESNORM = min(RESNORM_per_start_point);

% Calculate the number of times we have found the global minima
count = 0;
for n = 1 : starting_points
    if (RESNORM_per_start_point(n) - min_RESNORM <= accepted_deviation)
        count = count + 1;
    end
end

% Calculate the percentage of runs where the minimum RESNORM value 
% was found across all runs starting from different initial points.
percentage = count / starting_points;

% Determine the number of runs required to achieve a 95% confidence level 
% in finding the global minimum.
num_runs = log(0.05) / log(1-percentage);

end