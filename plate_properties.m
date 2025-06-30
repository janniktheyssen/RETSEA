function [p] = plate_properties(p, calc)
% Compute derived vibrational properties for two plates.
%
%   p = PLATE_PROPERTIES(p, calc) adds derived quantities to the input
%   structure array p for a two-plate system, based on material and 
%   geometric properties and the excitation frequency specified in calc.
%
%   Inputs:
%     p    - Struct array with material (E, E2, nu, eta), geometric (h, h2, A, perimeter),
%            and mass per unit area (m) properties for each plate.
%     calc - Struct with excitation frequency f (Hz).
%
%   Outputs:
%     p    - Updated struct array with additional fields: D, cb, cg, ma, kb, n,
%            lambda, geo.l, attenuation, count, k_dim, and Y_inf.
% 
% (c) Jannik Theyssen, LVA INSA Lyon, 2025 (GNU GPLv3)


for i = 1:2
    wc = calc.f * 2 * pi;
    % D - Bending stiffness of plate 1
    % ignoring the contribution of the damping
    % layer to the bending stiffness of plate 1
    p(i).D = p(i).mat.E * p(i).geo.h^3 / (12 * (1 - p(i).mat.nu^2)) + p(i).mat.E2 * p(i).geo.h2 * ((p(i).geo.h2+p(i).geo.h)/2)^2;
    % Bending wave speed (m/s)
    p(i).cb = (p(i).D / p(i).m)^(0.25) * sqrt(wc);
    % group speed (m/s)
    p(i).cg = 2 * p(i).cb;
    % attenuation coefficient (1/m)
    p(i).ma = p(i).mat.eta * wc / p(i).cg;
    % bending wavenumber (rad / m)
    p(i).kb = wc / p(i).cb;
    % modal density (1/Hz)
    p(i).n = p(i).geo.A / (4*pi) * sqrt(p(i).m/p(i).D);
    % wavelength (m)
    p(i).lambda = 2*pi / p(i).kb;
    % mean free path
    p(i).geo.l = pi * p(i).geo.A / p(i).geo.perimeter;
    % attenuation
    p(i).attenuation=p(i).mat.eta*wc./p(i).cg.*p(i).geo.l;
    % mode count
    p(i).count = p(i).n*wc; 
    % Dimensionless wavenumber
    p(i).k_dim = (p(i).m*wc^2./p(i).D).^0.25/2/pi.*p(i).geo.l;
    % Mobility of an infinite plate
    p(i).Y_inf = 1/(8 * sqrt(p(i).m * p(i).D));
    % Mean free path length per wavelength
    p(i).mfpl = p(i).geo.l / p(i).lambda;
end
