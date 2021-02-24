function [sub, varargout] = fun_io_octree_coordinate_to_array_sub(octree_coordinate)
% fun_io_octree_coordinate_to_array_sub converts the coordinate of the
% leaf in an 3D octree to its subscript in the cell array
% Input: 
%   octree_coordinate: N-by-1 numerical array. Each element is an integer
%   between 1 and 8, start from the first layer. 
% Output: 
%   sub: 3-by-1 integer array, subscript of the leaf in the cell array
%   ind: integer, index of the leaf in the cell array
% Note: 
%   This function is for handeling the rendered data from Janelia. The
%   data is organized in a row-major octree. 

    if isempty(octree_coordinate)
        sub = [1,1,1];
        if nargout > 1
            varargout{1} = 1;
        end
        return;
    end
    
    num_level = numel(octree_coordinate);
    cell_array_numel_line = 2^num_level;
    cell_array_numel_plane = cell_array_numel_line * cell_array_numel_line;
%     [ind_offset(:,1), ind_offset(:,2), ind_offset(:,3)] = ind2sub([2,2,2], 1:8);
%     ind_offset = tmp_offset - 2;
    % The following array is the same as the output of the two line above
    ind_offset = [-1,-1,-1;0,-1,-1;-1,0,-1;0,0,-1;-1,-1,0;0,-1,0;-1,0,0;0,0,0];
    % Pay attention to the order. 
    [tmp_j, tmp_i, tmp_k] = ind2sub([2,2,2], octree_coordinate(1));
    for idx = 2 : num_level
        tmp_i = tmp_i * 2 + ind_offset(octree_coordinate(idx),2);
        tmp_j = tmp_j * 2 + ind_offset(octree_coordinate(idx),1);
        tmp_k = tmp_k * 2 + ind_offset(octree_coordinate(idx),3);
    end
    sub = [tmp_i, tmp_j, tmp_k];
    if nargout > 1
        % Convert the subscrip to index
        varargout{1} = tmp_i + (tmp_j - 1) * cell_array_numel_line + (tmp_k - 1) * cell_array_numel_plane;
    end
end