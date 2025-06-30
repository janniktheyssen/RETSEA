function [Z] = receptance(plate, source_idx, receiver_idx, w)
%RECEPTANCE Compute receptance matrix via modal superposition.
%
%   Z = RECEPTANCE(plate, source_idx, receiver_idx, w) returns the 
%   receptance (frequency response function) between a set of source and 
%   receiver nodes on a vibrating plate structure using modal superposition.
%
%   Inputs:
%     plate        - Struct containing modal data with fields:
%                      .f0   - Vector of modal frequencies (Hz)
%                      .eta  - Modal damping factor
%                      .phi  - Mode shapes (num_modes x num_nodes)
%                      .norm - Modal normalization factors
%     source_idx   - Indices of source locations (scalar or vector)
%     receiver_idx - Indices of receiver locations (vector)
%     w            - Angular frequency (rad/s), scalar or vector
%
%   Output:
%     Z            - Receptance matrix of size [length(w) x length(receiver_idx)],
%                    giving the displacement at each receiver due to a unit force
%                    at source_idx, for each frequency in w.
% 
% (c) Jannik Theyssen, LVA INSA Lyon, 2025 (GNU GPLv3)

w0 = plate.f0 * 2 * pi; % eigenfrequencies of the modes (rad/s)

% precalculate the denominator for speed
denom = (plate.norm .* (w0.^2 + 1j*plate.eta.*w0.*w0 - w.^2));

Z = zeros(length(w), length(receiver_idx));

for i = 1:length(receiver_idx)
    Z(:,i) = sum(plate.phi(:,source_idx).*plate.phi(:,receiver_idx(i))./denom, 1);
end