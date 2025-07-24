%% ========================================================================
% AUTHOR:    Jannik Theyssen
% EMAIL:     jannik.theyssen@gmail.com
% AFFILIATION: Laboratory of Vibrations and Acoustics (LVA), INSA Lyon (FR)
% DATE:      2025-06-30
% VERSION:   v1.0
%
% DESCRIPTION:
%   Solves the hybrid radiative energy transfer - SEA model for two coupled
%   steel plates, one of  which is highly damped. 
%   Also solves the same problem using SEA and FEM via modal superposition. 
%   It then stores all three results for later comparison. 
%   This code is supplemental material to the article "A hybrid radiative 
%   transfer - SEA theory for point-coupled subsystems".
%
% LICENSE:
%   Copyright (c) 2025 Jannik Theyssen, LVA, INSA Lyon, FR.
% 
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY. See the GNU General Public License for 
%   details. The GNU General Public License is available at 
%   https://www.gnu.org/licenses/.
% 
% ------------------------------------------------------------------------
% For questions or feedback, please contact the author.
%% ========================================================================

clear; close all; clc

%% Geometries and materials
% geo_mat defines geometry and material properties of the two plates. These
% are stored in the "p" struct, which is loaded directly to the workspace.
geo_mat  

%% Calculation settings
% Calculation settings and results are saved in the "calc" struct array. 

% Enable plotting
calc(1).plotting = false;

% Reference model settings
calc(1).Nf = 6;  % number of frequency steps in half-power bandwidth

% Radiosity model settings
calc(1).ngp = 16;  % number of Gauss points
calc(1).element_size = 0.05;

% Spring stiffness in N/m
calc(1).K = 428;

% Center frequency in Hz
calc(1).f = 2000;

% Load field grid used in FreeFem to be used in ref and rad calculation
temp = load('eigenmodes_bunimovic_1mm.mat', 'coordinates');     
calc(1).field_grid = temp.coordinates(:, 1:end);
clear temp

% Define possible locations for the spring and source point. 
% The modal basis of the reference calculation is the limit for this resolution.
% Limit the selection of source and spring positions to a rectangle
% centered on the Bunimovic stadium
x_min = -0.04; 
x_max = 0.44;
y_min = p(1).geo.radius - 0.08;
y_max = p(1).geo.radius + 0.08;
point_grid_indices = find(...
    calc(1).field_grid(1,:) >= x_min & ...
    calc(1).field_grid(1,:) <= x_max & ...
    calc(1).field_grid(2,:) >= y_min & ...
    calc(1).field_grid(2,:) <= y_max);

eta = [0.001, 0.1];

% Generate all possible combinations of parameters
[spring_idx_grid, source_idx_grid, eta_grid] = ...
    ndgrid(point_grid_indices, point_grid_indices, eta);

% Convert set of parameters to a list
param_list = [spring_idx_grid(:), source_idx_grid(:), eta_grid(:)];

% Limit the parameter lsit to a random subset of all combinations
N_cases = 2200;
param_list = param_list(randperm(length(param_list), N_cases),:);

N_calc = size(param_list,1);

% Copy the general simulation settings to each case
calc = repmat(calc, [N_calc, 1]);


%% Run the simulations

disp('Running simulations...')
for ci = 1:N_calc
    tic
    % Add case-specific simulation settings
    calc(ci).source = calc(ci).field_grid(:, param_list(ci, 1))';
    calc(ci).spring_plate1 = calc(ci).field_grid(:, param_list(ci, 2))';

    % Choose a random connection point on plate 2
    calc(ci).spring_plate2 = calc(ci).field_grid(:, randsample(point_grid_indices,1))';

    disp(['Source   point on plate 1 at (', ...
        num2str(calc(ci).source, '%.2f '), ') m.'])
    disp(['Coupling point on plate 1 at (', ...
        num2str(calc(ci).spring_plate1, '%.2f '), ') m.'])
    disp(['Coupling point on plate 2 at (', ...
        num2str(calc(ci).spring_plate2, '%.2f '), ') m.'])

    p(1).mat.eta = param_list(ci,3);

    % Calculate the properties of each plate, store that in "p" struct, and
    % attach that to the calculation settings
    calc(ci).p = plate_properties(p, calc(ci));

    % Calculations
    calc(ci).ref = reference_model(calc(ci));
    calc(ci).sea = SEA(calc(ci));
    calc(ci).rad = radiative_energy_transfer(calc(ci));

    calc(ci).time = toc;
    
    fprintf('Case %d/%d done in %.2f s\n', ci, N_cases, calc(ci).time);
end

save('final_study.mat', 'calc')

%% Plot the results
postprocessing
