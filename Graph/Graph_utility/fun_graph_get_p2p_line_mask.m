function p2p_vec_ind = fun_graph_get_p2p_line_mask(mask_size, p1_sub, p2_sub)
% fun_graph_get_p2p_cylinder_mask construct a mask connecting two point.
% Input: 
%   mask_size: 3-by-1 numerical vector
%   p1_sub: 3-by-1 numerical vector, position of the point 1
%   p2_sub: 3-by-1 numerical vector, position of the point 2
%   r: radius of the mask 
% Output: 
%   mask: 3D logical array of size mask_size
% 

% Normalized direction vector from point 1 to point 2
num_dims = numel(p1_sub);
assert(num_dims == numel(mask_size), 'Dimension mismatch');
if iscolumn(p1_sub)
    p1_sub = p1_sub';
end
if iscolumn(p2_sub)
    p2_sub = p2_sub';
end

p2p_vec = p2_sub - p1_sub;
p2p_vec_norm = norm(p2p_vec);
p2p_vec_n = p2p_vec ./ p2p_vec_norm;

% p2p_vec_sub_floor = floor(p1_sub + p2p_vec_n .* (0 : 1 : p2p_vec_norm)');
% p2p_vec_sub_ceil = ceil(p1_sub + p2p_vec_n .* (0 : 1 : p2p_vec_norm)');
% p2p_vec_sub = cat(1, p2p_vec_sub_floor, p2p_vec_sub_ceil);
p2p_vec_sub = round(p1_sub + p2p_vec_n .* (0 : 1 : p2p_vec_norm)');
switch num_dims
    case 2
        p2p_vec_ind = sub2ind(mask_size, p2p_vec_sub(:, 1), p2p_vec_sub(:, 2));
    case 3
        p2p_vec_ind = sub2ind(mask_size, p2p_vec_sub(:, 1), p2p_vec_sub(:, 2), p2p_vec_sub(:, 3));
end
p2p_vec_ind = unique(p2p_vec_ind, 'stable');
end