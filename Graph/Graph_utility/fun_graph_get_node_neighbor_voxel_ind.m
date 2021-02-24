function neighbor_ind_list = fun_graph_get_node_neighbor_voxel_ind(node_ind, image_size, neighbor_type)
% fun_graph_get_node_neighbor_voxel_ind compute the linear indices for the
% voxels of the node and return it as a vector. 
% Input: 
%   node_ind: numerical vector, linear indices of the node voxel 
%   image_size: size of the 3D array. 
%   neighbor_type: 6 or 26 for 3D images
% Output: 
%   neighbor_ind_list: numerical vector
% 
% Implemented by Xiang Ji on 05/04/2019

if nargin < 3
    neighbor_type = 26;
end
assert(numel(image_size) == 3, 'The size of the image array should have 3 elements');

im_size_pad = image_size + 2;

persistent ind_adder image_size_used neighbor_type_used
if isempty(ind_adder) || any(image_size_used ~= im_size_pad) || neighbor_type_used ~= neighbor_type
    image_size_used  = im_size_pad;
    neighbor_type_used = neighbor_type;
    ind_adder = fun_skeleton_neighbor_add_coeff_3D(image_size_used, neighbor_type_used, true);
end
node_sub = fun_ind2sub(image_size, node_ind);
node_sub_pad = node_sub + 1;
node_ind_pad = sub2ind(im_size_pad, node_sub_pad(:, 1), node_sub_pad(:, 2), ...
    node_sub_pad(:, 3));
if iscolumn(node_ind_pad)
    node_ind_pad = node_ind_pad.';
end
if isrow(ind_adder)
    ind_adder = ind_adder.';
end
node_ind_pad_neighbor_ind = bsxfun(@plus, node_ind_pad, ind_adder);
node_neighbor_sub_pad = fun_ind2sub(im_size_pad, node_ind_pad_neighbor_ind(:));
node_neighbor_sub = node_neighbor_sub_pad - 1;
node_neighbor_valid_Q = all(bsxfun(@ge, node_neighbor_sub, 1), 2) & ...
    all(bsxfun(@le, node_neighbor_sub, image_size), 2);
% num_node_voxel = numel(node_ind);
neighbor_ind_list = sub2ind(image_size, node_neighbor_sub(node_neighbor_valid_Q, 1), node_neighbor_sub(node_neighbor_valid_Q, 2), node_neighbor_sub(node_neighbor_valid_Q, 3));
neighbor_ind_list = unique(neighbor_ind_list);
end

