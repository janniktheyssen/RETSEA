function [gamma] = view_factors_gamma(calc, elm)
%Compute view factors (gamma) for energy exchange between elements.
%
%   gamma = VIEW_FACTORS_GAMMA(calc, elm) calculates the view factor matrix 
%   for a set of elements, integrating over each linear element using Gaussian 
%   quadrature. The Gauss points act as discrete source points, and the element 
%   centers act as receiving points.
%   See equation (5.1) in accompanying publication.
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
%   Output:
%       gamma - View factor matrix [elm.n x elm.n], describing the 
%                  energy exchange factor between each source element (row)
%                  and receiving element (column).
%
%   The function performs the following steps:
%     1. Computes the distance (radius) between each Gauss point and each
%        element center.
%     2. Computes the cosine of the emission angle (theta) and incidence 
%        angle (phi) for each pair.
%     3. Evaluates the kernel G for each radius.
%     4. Integrates over Gauss points with Gaussian quadrature.
% 
% (c) Jannik Theyssen, LVA INSA Lyon, 2025 (GNU GPLv3)

% Number of Gauss points per element
n_gp = size(elm.gp, 2);

% Transpose and expand Gauss points to [n_source_elements x n_gp x 2]
gps = permute(elm.gp, [3, 2, 1]);

% Expand Gauss points [n_source_elements x n_gp x 1 x 2]
gps_exp = reshape(gps, [elm.n, n_gp, 1, 2]);

% Collocation points [1 x 1 x n_receiving_elements x 2]
centers = reshape(elm.centers, [1 1 elm.n 2]);

% Calculate vectors between source points and receiving element centers
% [n_source_elements x n_gp x n_cp x 2]
r_vec = gps_exp - centers;

% Squared distance and Euclidean norm [n_source_elements x n_gp x n_cp]
r_matrix = sqrt(sum(r_vec.^2, 4));

% Expand element normals to match source elements
elm_norm_source = reshape(elm.normals, [elm.n, 1, 1, 2]);

% Calculate cosine of the reflected angle via dot product
cos_theta = sum(elm_norm_source .* r_vec, 4) ./ r_matrix;

% Expand element normals to match receiver points
elm_norm_receiver = reshape(elm.normals, [1, 1, elm.n, 2]);

% Calculate incidence angle at the collocation points via dot product
cos_phi = sum(elm_norm_receiver .* (-r_vec), 4) ./ r_matrix;

% Evaluate the energy kernel function for all radii
G_matrix = G(r_matrix, calc);

% Gaussian integration by weighted sum
gamma = squeeze(sum(((G_matrix .* cos_theta .* cos_phi) .* elm.gw), 2)) / 2 .* repmat(elm.lengths,1,length(elm.elements));

end