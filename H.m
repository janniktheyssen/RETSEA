function H_phi = H(calc, elm, point)
% Compute the function H as defined in equation (2.2) in accompanying publication.
% 
%   H(calc, elm, point) calculates the energy density influence function at the 
%   centre of each boundary element due to the source located at point "point", 
%   considering the incidence angle phi at the boundary centre.
% 
%   Inputs: 
%       calc - Calculation struct containing, among others, geometry and material data
%
%       elm  - Structure describing the mesh elements, containing:
%                elm.n        - Number of elements
%                elm.gp       - Gauss points for each element (3D array: [2 x n_gp x elm.n])
%                elm.centers  - Center coordinates for each element [elm.n x 2]
%                elm.normals  - Inward normal vector for each element [elm.n x 2]
%                elm.gw       - Gauss weights (vector of length n_gp)
%                elm.lengths  - Length of each element (vector of length elm.n)
%                elm.elements - Element definition (nodes) [elm.nx2]
%
%       point - The source point coordinate [1x2]
%  
%   Outputs:
%       H_phi - The vector H [elm.n x 1]
% 
% (c) Jannik Theyssen, LVA INSA Lyon, 2025 (GNU GPLv3)

% Calculate the radius from the source point to all collocation points
r_vecs = elm.centers - point;
radius = sqrt(sum(abs(r_vecs).^2, 2));

% Calculate the cosine of the reflection angle by dot product of the 
% element normal vector and the vector from source point to element centre
cos_phi = sum(r_vecs .* elm.normals, 2) ./ radius;

% Calculate H vector
H_phi = G(radius, calc) .* cos_phi;

end