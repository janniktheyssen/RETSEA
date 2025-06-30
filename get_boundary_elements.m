function [elm, fig] = get_boundary_elements(elm, fig)
% Generate boundary elements and their properties from a closed 2D node list.
%   [elements, element_centers, element_normals, fig] = GET_BOUNDARY_ELEMENTS(nodes, fig)
%
%   Inputs:
%     nodes - Nx2 array
%         List of 2D coordinates defining the boundary (assumed to be ordered and closed).
%         Each row is a node [x, y]. The function will automatically close the loop.
%
%     fig - figure handle or false
%         If a figure handle is given (e.g., a number), the function will plot the boundary
%         elements, their centers, and inward normal vectors. If set to false, no plot is generated.
%
%   Outputs:
%     elements - Nx2 array
%         Each row defines a line element by specifying the indices of its two nodes.
%         The i-th row corresponds to element i, connecting nodes i and i+1 (looping back to 1 at the end).
%
%     element_centers - Nx2 array
%         Midpoints of each element, computed as the average of its two end nodes.
%
%     element_normals - Nx2 array
%         Inward-pointing unit normal vectors for each element, assuming the boundary is ordered
%         counter-clockwise. For each element, the normal is perpendicular to the element vector.
%
%     fig - figure handle or empty
%         If plotting is requested, the handle to the figure is returned.
%         Otherwise, returns [].
% 
% (c) Jannik Theyssen, LVA INSA Lyon, 2025 (GNU GPLv3)

    % Define elements and element center points
    n_nodes = size(elm.nodes, 1);
    elements = [(1:n_nodes)', [2:n_nodes, 1]']; % Define element connections
    element_centers = (elm.nodes(elements(:, 1), :) + elm.nodes(elements(:, 2), :)) / 2;

    % Calculate inward normal vectors
    element_normals = zeros(n_nodes, 2);
    for ei = 1:n_nodes
        dx = elm.nodes(elements(ei, 2), 1) - elm.nodes(elements(ei, 1), 1);
        dy = elm.nodes(elements(ei, 2), 2) - elm.nodes(elements(ei, 1), 2);
        normal = [dy, -dx] / norm([dy, -dx]);
        element_normals(ei, :) = normal;
    end

    % Plot if requested
    if fig ~= false
        figure(fig);
        hold on;
        % Plot nodes and elements
        % Plot element centres
        scatter(element_centers(:, 1), element_centers(:, 2), 'ro', 'filled', 'DisplayName', 'Element Centers');
        % Plot element normals
        quiver(element_centers(:, 1), element_centers(:, 2), ...
               element_normals(:, 1) * 0.002, element_normals(:, 2) * 0.002, ...
               'k', 'DisplayName', 'Element Normals');
        axis equal;
        legend;
        title('Boundary Elements and Normals');
        hold off;
        fig = gcf; % Return the figure handle
    else
        fig = [];
    end


    % Calculate the length of each element
    element_lengths = sqrt(sum((elm.nodes(elements(:, 2), :) - elm.nodes(elements(:, 1), :)).^2, 2));


    elm.elements = elements;
    elm.centers = element_centers;
    elm.normals = element_normals;
    elm.lengths = element_lengths;
    elm.n = length(elements);

end
