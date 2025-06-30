function W_fp = field_evaluation(calc, elm, sigma)
%Compute the energy density W_fp in the domain due to boundary sources.
% 
%  Inputs:
%    calc: struct with problem data (contains field_grid, p, etc.)
%    elm:  struct with element data (gp, normals, lengths, gw, etc.)
%    sigma: boundary source strengths for each element
%  Output:
%    W_fp: [n_fp x 1] energy density at each field point
% 
% 
% (c) Jannik Theyssen, LVA INSA Lyon, 2025 (GNU GPLv3)

field_grid = calc.field_grid';

n_fp = size(field_grid, 1);
n_S = elm.n;
n_gp = calc.ngp;

% Expand Gauss points to [n_S x n_gp x 1 x 2]
gp = permute(elm.gp, [3, 2, 1]);
gp = reshape(gp, [n_S, n_gp, 1, 2]);

% Expand field points [1 x 1 x n_fp x 2]
fp = reshape(field_grid, [1, 1, n_fp, 2]);

% Vector from Gauss points to field points [n_S x n_gp x n_fp x 2]
vec = gp - fp;

% Magnitude of vectors [n_S x n_gp x n_fp]
radius = sqrt(sum(vec.^2, 4)); 

% Normalise vectors [n_S x n_gp x n_fp x 2]
norm_vec = vec ./ radius; 

% Expand normals (normal on each source element) [n_S x 1 x 1 x 2]
normals = reshape(elm.normals, [n_S, 1, 1, 2]);

% Dot product cos(phi) [n_S x n_gp x n_fp]
cos_phi = sum(normals .* norm_vec, 4);

% Evaluate G for each radius [n_S x n_gp x n_fp]
G_vals = G(radius, calc);

% Integrate: sum over Gauss points for each element → then sum over elements
% Gauss weights: [1 x n_gp]
gw = reshape(elm.gw, [1, n_gp]);

% Contribution from each boundary element to each field point [n_S x n_fp]
W_fp = squeeze(sum( G_vals .* cos_phi .* gw, 2) ./ 2 ) .* sigma .* elm.lengths; 

% Sum over elements, transpose to [n_fp x 1]
W_fp = sum(W_fp, 1)';
end