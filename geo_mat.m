% GEO_MAT define relevant geometry and material properties for the
% simulation.
% 
% (c) Jannik Theyssen, LVA INSA Lyon, 2025 (GNU GPLv3)

%% Geometry definition
% The geometry of the Bunimovic stadium is identical to the one published
% at https://royalsocietypublishing.org/doi/10.1098/rspa.2021.0488, Fig. 1b)

% COORDINATES
% Define the stadium geometry by its radius and length. 
% 
%                        _
%     xxxxxxxxxxxxxx     ^
%   xx              xx   | radius
%  xx                xx  |  
% x                    x v
% x   ^                x 
% x   | y              x
% x   |                x
%  xx |              xx 
%   xx|             xx  
%     .--->x xxxxxxx
%    (0,0)
%     |<- length ->|
% |<----- width  ----->|

geo.diameter = 0.378; % Diameter (m)
geo.radius = geo.diameter/2;   % Radius, see sketch below. 
geo.length = 0.4;  % length of central rectangle  (m)
geo.width = geo.length+2*geo.radius; 
% surface area of the stadium
geo.A = (geo.width-2*geo.radius) * 2*geo.radius + geo.radius^2 * pi;   
geo.perimeter = 2*geo.length + geo.diameter*pi;

% Both plates share the surface geometry
p(1).geo = geo;
p(2).geo = geo;
clear geo;

% h is the plate thickness (m)
p(1).geo.h = 0.002;
p(2).geo.h = 0.001;

% h2 is the thickness of the damping layer (m)
p(1).geo.h2 = 0.004;
p(2).geo.h2 = 0;


%% Material properties

steel.E = 210e9;   % Young's modulus of the carrier plate (N/m^2)
steel.rho = 7800;  % Density of steel (kg/m^3)
steel.nu = 0.3;    % Poisson ratio

p(1).mat = steel;
p(1).mat.E2 = 210e6; % Young's modulus of the damping layer plate 1 (N/m^2)
p(2).mat = steel;
p(2).mat.E2 = 0; % Young's modulus of the damping layer plate 2 (N/m^2)
clear steel

p(1).mat.eta = 0.1;    % Loss factor in plate 1    
p(2).mat.eta = 0.0005; % Loss factor in plate 2

p(1).m = p(1).geo.h * p(1).mat.rho + 4.5; % Mass per unit area of plate 1 (kg/m^2)
p(2).m = p(2).geo.h * p(2).mat.rho;       % Mass per unit area of plate 2 (kg/m^2)

p(1).M = p(1).m * p(1).geo.A;         % Total mass of plate 1
p(2).M = p(2).m * p(2).geo.A;         % Total mass of plate 2

% Reference model calculation
p(1).freefem_result_file = 'eigenmodes_bunimovic_2mm_damped.mat';
p(2).freefem_result_file = 'eigenmodes_bunimovic_1mm.mat';
