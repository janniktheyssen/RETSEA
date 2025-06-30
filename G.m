function [G] = G(r, calc)
% Compute the energy kernel (Green's function-like quantity) 
% based on input distance vector r
%
%   G = G(r, calc) computes the function G for each element in the vector r.
%   The computation follows a piecewise definition (see Eq. (4.9)),
%   applying different formulas depending on the value of r relative to the 
%   wavelength lambda.
% 
% (c) Jannik Theyssen, LVA INSA Lyon, 2025 (GNU GPLv3)

% Initialize G output vector
G = zeros(size(r));

% Determine the indices for each condition
idx1 = r > calc.p(1).lambda / pi^2;
idx2 = ~idx1;

% Define the two functions
f1 = @(r) 1 ./ (2*pi*calc.p(1).cg*r) .* exp(-calc.p(1).ma*r);
f2 = @(r) sqrt(calc.p(1).mat.rho*calc.p(1).geo.h/calc.p(1).D) / 8;

% Apply the functions selectively
G(idx1) = f1(r(idx1));
G(idx2) = f2(r(idx2));