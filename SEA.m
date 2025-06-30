function res = SEA(calc)
% Estimate energy distribution using SEA between two plates.
%   res = SEA(calc) computes the vibrational energy on two 
%   coupled plates using SEA based on modal densities, damping, and coupling stiffness.
%
%   Inputs:
%     calc - Struct with fields:
%            f - Center frequency of the octave band (Hz)
%            K - Coupling stiffness between the plates (N/m)
%            p - Struct array of length 2, each with fields:
%                m - Mass of the plate (kg/m2)
%                D - Flexural rigidity of the plate
%                n - Modal density (modes per rad/s)
%                mat.eta - Loss factor (dimensionless)
%
%   Output:
%     res.energy - 2x1 vector of estimated vibrational energy on each plate
%     res.beta   - Coupling loss factor between the plates
%
%   Reference:
%     For the coupling loss factor beta, see
%     B.R. Mace and L. Ji, "The statistical energy analysis of coupled sets 
%     of oscillators", Proc. R. Soc. A, 463 (2007), 1359–1377.
% 
% (c) Jannik Theyssen, LVA INSA Lyon, 2025 (GNU GPLv3)

p = calc.p;

wc = calc.f * 2 * pi;

% coupling power proportionality factor
betaSEA = calc.K^2/(32 * pi * wc^2 * sqrt(p(1).m * p(1).D * p(2).m * p(2).D));

SEA_A_matrix = [p(1).mat.eta*wc + betaSEA/p(1).n, -betaSEA/p(2).n;
     -betaSEA/p(1).n, p(2).mat.eta*wc + betaSEA/p(2).n];

Win = 1; % Unit power input to plate 1

res.energy = SEA_A_matrix\[Win;0];
res.beta = betaSEA;
