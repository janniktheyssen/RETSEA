function res = reference_model(calc)
% Compute reference vibro-acoustic response using modal data from FreeFEM++
%
%   res = REFERENCE_MODEL(p, calc) calculates the energy density field on 
%   two plates and power transfer between these two plates in a band-
%   limited frequency range. The calculation is uses a modal 
%   superposition of a FreeFEM++ modal analysis result.
%
%   Inputs:
%     p    - Struct array (length 2) with fields:
%              .freefem_result_file: path to FreeFem modal result file
%              .mat: material properties (including eta)
%              .geo: geometry (including area A)
%              .m: mass per unit area
%     calc - Struct with calculation parameters:
%              .f: central frequency (Hz)
%              .Nf: number of frequency steps
%              .K: spring stiffness between plates
%              .source: [x, y] coordinates of excitation
%              .spring_plate1/2: [x, y] coordinates of coupling points
%              .field_grid: [2 x N] grid of receiver positions
%              .plotting: flag for coordinate index visualization
%
%   Outputs:
%     res - Struct with fields:
%              .energy_field: energy density at field points [2 x N]
%              .energy: total energy on each plate
%              .beta: energy transmission coefficient
%              .power: total transmitted power
%              .Y_in: input mobility over frequency
%              .P_in: input power
%              .P_diss: dissipated power on each plate
% 
% (c) Jannik Theyssen, LVA INSA Lyon, 2025 (GNU GPLv3)

p = calc.p;
% Load FreeFem calculation
for i = 1:2
    p(i).ff = load(p(i).freefem_result_file, 'f0', 'norm', 'phi', 'coordinates');
    p(i).ff.w0 = 2 * pi * p(i).ff.f0;
    p(i).ff.eta = p(i).mat.eta;
end

% Plot grid
if calc.plotting
    figure(101)
    text(p(1).ff.coordinates(1,:), p(1).ff.coordinates(2,:), num2cell(1:length(p(1).ff.coordinates)))
    xlim([-0.25, 0.65])
    ylim([-0.1, max(p(1).ff.coordinates(2,:))+0.1])
    title('Field coordinate indices in reference calculation')
end

wc = 2*pi*calc.f;             % Angular frequency (rad/s)
wmin=wc/sqrt(2);              % min frequency of octave
wmax=wc*sqrt(2);              % max frequency of octave

mat = [p.mat];                % extract material data from plate struct

% Calculate exact modal density
for i = 1:2
    p(i).n_exact = sum(p(i).ff.w0 < wmax & p(i).ff.w0 > wmin)/(wmax-wmin);
end
dw=min(min(1./[p.n_exact],[mat.eta]*wc))/calc.Nf;  % frequency step
ww=wmin:dw:wmax;      % frequency vector

% Limit FreeFem modes up to N*w_max
N = 2;
for i = 1:2
    widx = p(i).ff.w0 < (wmax * N);
    p(i).ff.w0  = p(i).ff.w0(widx);
    p(i).ff.norm  = p(i).ff.norm(widx);
    p(i).ff.phi = p(i).ff.phi(widx,:);
    p(i).ff.f0 = p(i).ff.f0(widx);
end

% Definition of source and connection point coordinates
r1s_idx = find_grid_point_idx(p(1).ff.coordinates, calc.source(1), calc.source(2));
r1c_idx = find_grid_point_idx(p(1).ff.coordinates, calc.spring_plate1(1), calc.spring_plate1(2));
r2c_idx = find_grid_point_idx(p(1).ff.coordinates, calc.spring_plate2(1), calc.spring_plate2(2));

r1r_idx = find_grid_point_idx(p(1).ff.coordinates, calc.field_grid(1,:), calc.field_grid(2,:));
r2r_idx = find_grid_point_idx(p(2).ff.coordinates, calc.field_grid(1,:), calc.field_grid(2,:));


% Precalculate relevant receptances to make modal superpostion more
% efficient
r1c = receptance(p(1).ff, r1c_idx, r1c_idx, ww);   % input receptance at spring position plate 1
r2c = receptance(p(2).ff, r2c_idx, r2c_idx, ww);   % input receptance at spring position plate 2
r1sc = receptance(p(1).ff, r1s_idx, r1c_idx, ww);  % source to spring plate 1
r1sr = receptance(p(1).ff, r1s_idx, r1r_idx, ww);  % source to receiver plate 1
r1cr = receptance(p(1).ff, r1c_idx, r1r_idx, ww);  % spring to reciver plate 1
r2cr = receptance(p(2).ff, r2c_idx, r2r_idx, ww);  % spring to receiver plate 2

% Preallocate field energy density vector
energy_field = zeros(2, length(r2r_idx));

% Initialise total power
power = 0;

% Preallocate input mobility vector
Y_in = zeros(length(ww),1);

for iw = 1:length(ww)
    w=ww(iw);
    phi = [r1c(iw), -r1c(iw); 
           -r2c(iw),  r2c(iw)];

    Y_in(iw) = 1j*w*receptance(p(1).ff, r1s_idx, r1s_idx, w);
    F = 1;

    F0 = [F * r1sc(iw); 0];

    u = (eye(2)+calc.K*phi)\F0;

    H1r = r1sr(iw,:) + calc.K * (u(2)-u(1)) * r1cr(iw,:);
    H2r = - calc.K * (u(2)-u(1)) * r2cr(iw,:);

    H = [H1r;H2r];
    
    P_in = (wmax - wmin) / pi * p(1).Y_inf;
    % estimating the kinetic energy density (twice, assuming potential energy = kinetic energy)
    energy_field = energy_field + abs(H.^2) * w^2 .* [p.m]' * dw / pi / P_in;
    power = power + calc.K / pi * real(1i*w*conj(u(1))*u(2)) * dw / P_in;
end


% store results in res struct
res.energy_field = energy_field;  % energy density in the field
res.energy = mean(res.energy_field,2) * p(1).geo.A;  % total energy per plate (integration of the energy density over the surface A)
res.beta = power ./ (res.energy(1)/p(1).n - res.energy(2)/p(2).n);
res.power = power;
res.Y_in = Y_in;
res.P_in =  P_in;
res.P_diss(1) = (p(1).mat.eta * wc * mean(energy_field(1,:)) * p(1).geo.A);
res.P_diss(2) = (p(2).mat.eta * wc * mean(energy_field(2,:)) * p(2).geo.A) / res.power;

end
