function idx = find_grid_point_idx(grid, x_coord, y_coord)
% Find indices of points in a 2D grid matching given coordinates.
%
%   idx = FIND_GRID_POINT_IDX(grid, x_coord, y_coord) returns the indices of
%   grid points that match the given x and y coordinates within a tolerance.
%   The grid is assumed to be a 2D array of point coordinates with size [N x 2]
%   or [2 x N] and will be transposed if necessary.
%
%   Inputs:
%     grid     - Nx2 or 2xN array of [x, y] coordinates
%     x_coord  - Vector of x-coordinates to locate
%     y_coord  - Vector of y-coordinates to locate
%
%   Output:
%     idx      - Vector of indices in the grid corresponding to the input
%                (x_coord, y_coord) pairs
%
%   Notes:
%     - A fixed tolerance of 1e-4 is used to match coordinates.
%     - Each (x, y) pair must uniquely match a grid point.
% 
% (c) Jannik Theyssen, LVA INSA Lyon, 2025 (GNU GPLv3)


if size(grid,2) > size(grid,1)  % check orientation 
    grid = grid';
end
if size(x_coord,2) > size(x_coord,1)  % check orientation 
    x_coord = x_coord';
end
if size(y_coord,2) > size(y_coord,1)  % check orientation 
    y_coord = y_coord';
end
if size(grid,2) > 2
    error('Grid has more than two dimensions')
end

idx = zeros(length(x_coord), 1);

for i = 1:length(x_coord) 
    gridpoint_idx_x = find(abs(grid(:,1)-x_coord(i)) < 0.0001);
    grid_idx_y = find(abs(grid(gridpoint_idx_x,2) - y_coord(i)) < 0.0001);
    idx(i) = gridpoint_idx_x(grid_idx_y);
end

