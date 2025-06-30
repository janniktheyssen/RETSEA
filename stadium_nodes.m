function [elm, fig] = stadium_nodes(radius, length, distance, plot_flag)
    % Define the stadium geometry by its radius and length.
    % Returns a number of points spaced along the edges of the stadium 
    % spaced by `distance`.
    % 
    % (c) Jannik Theyssen, LVA INSA Lyon, 2025 (GNU GPLv3)

    % Length of the straight section
    l_straight_section = (length - 2 * radius);
    % Total circumference
    l_boundary = 2 * pi * radius + 2 * l_straight_section;

    % Number of nodes on the half-circle
    n_nodes_on_half_circle = ceil(pi * radius / distance);
    % Number of nodes on the straight section
    n_nodes_on_straight_section = round(l_straight_section / distance);
    % Total number of nodes
    n_nodes = (n_nodes_on_half_circle + n_nodes_on_straight_section) * 2;

    % Angles to the points in the circle
    phi = (0:(n_nodes_on_half_circle-1)) / n_nodes_on_half_circle * pi;

    % Position of the points in the straight sections
    x_straight = (0:(n_nodes_on_straight_section-1)) / n_nodes_on_straight_section * l_straight_section;

    % Nodes on the left circle
    nodes_left_circle = [radius * sin(-phi); radius * sin(phi+pi/2)+radius];
    % Nodes on the lower straight section
    nodes_lower_straight = [x_straight; -radius * ones(1, n_nodes_on_straight_section)+radius];
    % Nodes on the right circle
    nodes_right_circle = [radius * sin(phi) + l_straight_section; radius * sin(phi-pi/2)+radius];
    % Nodes on the upper straight section
    nodes_upper_straight = [l_straight_section - x_straight; radius * ones(1, n_nodes_on_straight_section)+radius];

    % Combine all nodes
    elm.nodes = [nodes_left_circle, nodes_lower_straight, nodes_right_circle, nodes_upper_straight]';

    % fprintf("Generated geometry with %d nodes, each %.2f cm apart.\n", n_nodes, distance);

    % Plot the geometry if requested
    if plot_flag
        figure;
        hold on;
        plot(nodes_left_circle(1,:), nodes_left_circle(2,:), 'x-', 'DisplayName', 'Left Circle');
        plot(nodes_lower_straight(1,:), nodes_lower_straight(2,:), 'x-', 'DisplayName', 'Lower Straight');
        plot(nodes_right_circle(1,:), nodes_right_circle(2,:), 'x-', 'DisplayName', 'Right Circle');
        plot(nodes_upper_straight(1,:), nodes_upper_straight(2,:), 'x-', 'DisplayName', 'Upper Straight');
        axis equal;
        legend;
        title('Stadium Geometry');
        hold off;
    end

    % Return an empty figure if plot_flag is false
    if ~plot_flag
        fig = false;
    else
        fig = gcf;
    end
end
