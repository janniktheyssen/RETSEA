function res = radiative_energy_transfer(calc)
% RADIATIVE_ENERGY_TRANSFER 
% res = radiative_energy_transfer(calc) Computes the energy 
% transfer by between two coupled plates where the first plate is a
% radiative subsystem and the second plate is an SEA subsystem
%
% INPUT:
%   calc - struct with fields:
%          .f (Hz): center frequency
%          .p: struct array with plate properties
%          .element_size (m): desired length for linear boundary elements
%          .plotting (bool): enable/disable plots
%          .ngp: number of Gauss points (currently only 16 supported)
%          .source: [x,y] source position
%          .spring_plate1: [x,y] spring position on plate 1
%          .spring_plate2: [x,y] spring position on plate 2 (not used directly here?)
%          .K (N/m): spring stiffness
%          .field_grid: [2 x N] grid of field points
%
% OUTPUT:
%   res - struct with fields:
%         .p_diss: dissipated power in plate 1 (W)
%         .energy: [E1; E2] energies in subsystem 1 and 2
%         .power: coupling power P12
%         .energy_field: energy density at each field point on plate 1
%         .beta: effective coupling factor
% 
% (c) Jannik Theyssen, LVA INSA Lyon, 2025 (GNU GPLv3)

% Centre frequency wc
wc = calc.f * 2 * pi; 

% Generate geometry
[elm, fig] = stadium_nodes(calc.p(1).geo.radius, calc.p(1).geo.width, calc.element_size, calc.plotting);

% Generate boundary elements, element centers, and element normals
[elm, fig] = get_boundary_elements(elm, fig);

% Calculate Gauss points and weights (only 16 GPs implemented)
[elm, ~] = get_gauss_points_and_weights(elm, fig, calc.ngp);

if calc.plotting
    figure(11);clf
    for i = 1:length(elm.nodes)
        plot([min(elm.nodes(:,1)), max(elm.nodes(:,1))], [0, 0])
        hold on
        text(elm.nodes(i,1), elm.nodes(i,2), num2str(i))
    end
    axis equal
end

% Set up system of equations (equation 5.4 in accompanying publication):
% 1. evaluate view factor matrix
gamma = view_factors_gamma(calc, elm);

% 2. evaluate Psi_phi(p_n, c) vector
Psi_phi_pc = H(calc, elm, calc.spring_plate1);

% 3. evaluate xi_n vector
xi_spring = xi(calc, elm, calc.spring_plate1);

% 4. evaluate right hand side Psi_phi(p_n, s)
Psi_phi_ps = H(calc, elm, calc.source);

% Calculate SEA coupling proportionality factor beta
betaSEA = calc.K^2/(32 * pi * wc^2 * sqrt(calc.p(1).m * calc.p(1).D * calc.p(2).m * calc.p(2).D));
% See B.R. Mace, L. Ji, The statistical energy analysis of coupled sets of oscillators, 
% Proceedings of the royal society A 463 (2007) 1359-1377.

% Some helper variables
bA1n1 = betaSEA*calc.p(1).geo.A/calc.p(1).n;
bn2 = betaSEA/calc.p(2).n;

% Assembling the matrix
A = [...
    % Radiosity power balance on the boundary (eq. 5.1)
    eye(elm.n)/(pi*calc.p(1).cg) - gamma, bA1n1*Psi_phi_pc, -bn2*Psi_phi_pc; ...
    % Power densiity in coupling point (eq. 5.2)
    -xi_spring, 1+bA1n1*G(0,calc), -bn2*G(0,calc); ...
    % SEA power balance (eq. 5.3)
    zeros(1, elm.n), -bA1n1, bn2+calc.p(2).mat.eta*wc];

% Assembling the right hand side
r_cs = sqrt(sum(abs(calc.source - calc.spring_plate1).^2));
P_in = 1;
b = P_in * [Psi_phi_ps; G(r_cs, calc); 0];

% Solving the matrix system
x = A\b;

% Extracting values from vector of unknowns
sigma = x(1:end-2); 
W1 = x(end-1);
E2 = x(end);

% Calculate coupling power P12 according to equation (3.1)
P12 = betaSEA * (W1 * calc.p(1).geo.A/calc.p(1).n - E2 / calc.p(2).n);

% Postprocessing: Calculate energy density at field points 
field_grid = calc.field_grid'; % [n_fp x 2]

% Direct contribution from source
field_vectors_source = calc.source - field_grid;  % [n_fp x 2]
field_radius_source = sqrt(sum(field_vectors_source.^2, 2));  % [n_fp x 1]
W_direct = G(field_radius_source, calc) * P_in;   % [n_fp x 1]

% Direct contribution from spring
field_vectors_spring = calc.spring_plate1 - field_grid;  % [n_fp x 2]
field_radius_spring = sqrt(sum(field_vectors_spring.^2, 2)); % [n_fp x 1]
% negative contribution as P12 is positive for loss in subsystem 1!
W_direct_spring = - G(field_radius_spring, calc) * P12;  % [n_fp x 1]


% Calculate contribution from edge reflections
W_fp = field_evaluation(calc, elm, sigma);

% Total field energy density
W = W_fp + W_direct + W_direct_spring;

% Store in result struct
res.p_diss = calc.p(1).mat.eta * wc * mean(W) * calc.p(1).geo.A;
res.energy = [mean(W)*calc.p(1).geo.A; E2];
res.power = P12;
res.energy_field = W;
res.beta = res.power./(res.energy(1)/calc.p(1).n-res.energy(2)/calc.p(2).n);
res.Wc = W1;