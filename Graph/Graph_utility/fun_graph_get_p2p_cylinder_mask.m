function mask = fun_graph_get_p2p_cylinder_mask(mask_size, p1_sub, p2_sub, r)
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
if iscolumn(p1_sub)
    p1_sub = p1_sub';
end
if iscolumn(p2_sub)
    p2_sub = p2_sub';
end

p2p_vec = p2_sub - p1_sub;
p2p_vec = p2p_vec ./ norm(p2p_vec);

[sub_1, sub_2, sub_3 ] = ndgrid(1:mask_size(1), 1:mask_size(2), 1:mask_size(3));
% Vector: position in the mask - point 2
vec_r_ep2 = cat(2, sub_1(:) - p2_sub(1), sub_2(:) - p2_sub(2), ...
    sub_3(:) - p2_sub(3));
vec_r_ep1 = cat(2, sub_1(:) - p1_sub(1), sub_2(:) - p1_sub(2), ...
    sub_3(:) - p1_sub(3));


% The distance^2 between the line connecting point 1 and point 2 to the
% true voxel point should be small than r^2
dp_r2ep2_p2p = vec_r_ep2 * p2p_vec';
tmp_dist = sum(vec_r_ep2.^2, 2) - (dp_r2ep2_p2p) .^ 2;
dp_r2ep1_p2p = vec_r_ep1 * p2p_vec';
mask = (tmp_dist <= (r*r)) & (dp_r2ep1_p2p >= 0) & (dp_r2ep2_p2p <= 0);
mask = reshape(mask, mask_size);
end