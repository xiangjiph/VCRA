function [ind, varargout] = fun_io_octree_coordinate_to_child_pos(octree_coordinate, child_level)
% fun_io_octree_coordinate_to_child_pos
% Input: 
%   octree_coordinate: numerical vector, specify the position of the leaf
%   in the octree
%   child_level: integer, specifying the level at which the child indices
%   and subscripts of the input leaf are calculated
% Output: 
%   ind: indices of all the sub-blocks at $child_level$ of the given leaf
%   varargout: subscript of these sub-blocks

parent_cell_sub = fun_io_octree_coordinate_to_array_sub(octree_coordinate);
parent_octree_level = numel(octree_coordinate);
child_cell_array_size = 2^(child_level - parent_octree_level);
[child_cell_sub_1, child_cell_sub_2, child_cell_sub_3] = ndgrid(child_cell_array_size* (parent_cell_sub(1) - 1) + 1 : child_cell_array_size * parent_cell_sub(1), ...
    child_cell_array_size* (parent_cell_sub(2) - 1) + 1 : child_cell_array_size * parent_cell_sub(2), ...
    child_cell_array_size* (parent_cell_sub(3) - 1) + 1 : child_cell_array_size * parent_cell_sub(3));
ind = sub2ind(ones(3,1)*(2^child_level), child_cell_sub_1(:), child_cell_sub_2(:), child_cell_sub_3(:));
if nargout > 1
    varargout{1} = [child_cell_sub_1(:), child_cell_sub_2(:), child_cell_sub_3(:)];
end


end