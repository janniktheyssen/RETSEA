function[xi_vec] = xi(calc, elm, point)
%Compute the xi function for a given source point (equation 5.2 in
%   accompanying publication).
% 
%   xi_vec = XI(calc, elm, point) integrates the energy kernel
%   times the cosine of the reflection angle between the inward normal and the 
%   vector between the element centre and the point "point"
%   over all elements and Gauss points.
%
%   Inputs:
%     calc  - Structure with parameters used by G()
%     elm   - Structure with:
%               elm.n        - Number of elements
%               elm.ngp      - Number of Gauss points per element
%               elm.gp       - Gauss points [2 x ngp x n]
%               elm.normals  - Normals [n x 2]
%               elm.gw       - Gauss weights [ngp x 1]
%               elm.lengths  - Element lengths [n x 1]
%     point - Receiver point [1 x 2]
%
%   Output:
%     xi_vec - The integrated xi value
% 
% (c) Jannik Theyssen, LVA INSA Lyon, 2025 (GNU GPLv3)


% Vector from source point to each Gauss point [n x ngp x 2]
r_vecs = permute(elm.gp, [3, 2, 1]) - reshape(point, 1, 1, 2);

% CAlculate the radius [n x ngp]
radius = sqrt(sum(r_vecs.^2, 3)); 

% Normalise vectors
r_mags = radius;
norm_vecs = r_vecs ./ r_mags;

% Calculate cosine of reflection angle theta
% Expand normals for broadcasting: 
normals_exp = reshape(elm.normals, elm.n, 1, 2);
cos_theta = sum(normals_exp .* norm_vecs, 3); % [n x ngp]

% Evaluate energy kernel function [n x ngp]
G_val = G(radius, calc);

% Integrate with Gauss weights and element lengths [1 x n]
xi_vec = ((G_val .* cos_theta) * elm.gw' / 2 .* elm.lengths)';

end
