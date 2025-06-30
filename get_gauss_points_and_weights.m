function [elm, fig] = get_gauss_points_and_weights(elm, fig, ngp)
% Compute Gauss points and corresponding weights for each element. 
% 
% [elm, fig] = get_gauss_points_and_weights(elm, fig, ngp) returns the
% Gauss points in global coordinates
%   Input:  elm - element struct
%           fig - figure handle
%           ngp - number of Gauss points (integer) NOTE: At the moment,
%           only 16 Gauss points are implemented.
%   Output:
%           elm - element struct, now including the 
%                       - Gauss points in global coordinates as elm.gp
%                       - Weights as elm.gw
%           fig - the figure handle, now including Gauss points
% 
% (c) Jannik Theyssen, LVA INSA Lyon, 2025 (GNU GPLv3)

% Evaluate Gauss points in local coordinates
[xi_gp, elm.gw] = Gauss(ngp);

gauss_points = zeros(2, ngp, elm.n);

% Compute Gauss point positions for each element (global coordinates)
for ei = 1:elm.n
    a = elm.nodes(elm.elements(ei, 1), :);
    b = elm.nodes(elm.elements(ei, 2), :);
    for gpi = 1:ngp
        gauss_points(:, gpi, ei) = ((b - a) / 2) * xi_gp(gpi) + (a + b) / 2;
    end
end

% Store global positions of Gauss points in elm struct
elm.gp = gauss_points;

% Plot Gauss points if requested
if fig ~= false
    figure(fig);
    hold on;
    % Reshape Gauss points for plotting
    reshaped_gp = permute(gauss_points, [3,2,1]);
    scatter(reshaped_gp(:,:, 1), reshaped_gp(:, :, 2), 'DisplayName', 'Gauss Points');
    % plot(nodes(:, 1), nodes(:, 2), 'bo-', 'DisplayName', 'Boundary Nodes');
    axis equal;
    legend;
    title('Gauss Points');
    hold off;
    fig = gcf; % Return the figure handle
else
    fig = [];
end

end

function [x, w] = Gauss(n)
%Gauss(n)
% 
% Return the Gauss points and weights for Gaussian quadrature
% 
%   Input: 
%       n: number of Gauss points (currently only 16 implemented)
%   Output: 
%       x: Abscissa
%       w: Weight

 switch n
     case 16
     x = [-0.989400934991650
    -0.944575023073233
    -0.865631202387832
    -0.755404408355003
    -0.617876244402644
    -0.458016777657227
    -0.281603550779259
    -0.0950125098376374
    0.0950125098376374
    0.281603550779259
    0.458016777657227
    0.617876244402644
    0.755404408355003
    0.865631202387832
    0.944575023073233
    0.989400934991650];
    
    w = [
    0.0271524594117541
    0.0622535239386478
    0.0951585116824929
    0.124628971255534
    0.149595988816577
    0.169156519395003
    0.182603415044924
    0.189450610455069
    0.189450610455069
    0.182603415044924
    0.169156519395003
    0.149595988816577
    0.124628971255534
    0.0951585116824929
    0.0622535239386478
    0.0271524594117541]';

     otherwise
         error('Number of Gauss integration points not implemented.')
 end
end