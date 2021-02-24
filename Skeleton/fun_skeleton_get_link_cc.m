function [link_cc, varargout] = fun_skeleton_get_link_cc(voxel_list, voxel_neighbor_idx_list, ...
    l_link_start_voxel_idx, voxel_unvisited_Q)
% fun_skeleton_get_link_cc is a subfunction of fun_skeleton_to_graph, which
% finds the link connected components. Except for the beginning voxels,
% each voxels in the link has exactly two neigbhors, whose indices in the
% voxel list in the first two row of voxel_neighbor_idx_list. 
% Input: 
%   voxel_list: linear indices of the voxels in the 3D space
%   voxel_neighbor_idx_list: 26-by-N array, where N is the number of voxels
%   l_link_start_voxel_idx: list index of the voxels to start tracking the
%   link connected component
%   voxel_unvisited: N-by-1 logical vector, where 1 means the voxel is not
%   a isolated point nor a node voxel
% Output: 
%   link_cc: structure with fields PixelIdxList, NumObjects and
%   Connectivity
%   num_unvisited_points: number of remianing unvisited voxels
%   voxel_unvisited_Q: N-by-1 logical vecotors, updated vector
% Implemented by Xiang Ji on 08/13/2019

[num_neighbors, num_voxel] = size(voxel_neighbor_idx_list);
assert(num_voxel == numel(voxel_list), 'Number of voxels in the voxel list is not the same as in the voxel_neighbor_idx_list');
num_unvisited_points = nnz(voxel_unvisited_Q);

link_cc = struct;
link_cc.PixelIdxList = cell(num_unvisited_points, 1);
link_cc.NumObjects = 0;
link_cc.Connectivity = num_neighbors;
if num_unvisited_points == 0
    if nargout > 1
        if nargout == 2
            varargout{1} = num_unvisited_points;
        elseif nargout == 3
            varargout{1} = num_unvisited_points;
            varargout{2} = voxel_unvisited_Q;
        end
    end
    return;
end
t_num_link_start_voxel = numel(l_link_start_voxel_idx);
t_link_idx = zeros(2000,1); % 2000um is a reasonable maximum length of the vessel segement you could findin the mouse brian. 
t_start_search_idx = 1;
t_num_cc = 0;
t_start_point_idx = 0;
keep_tracking = false;
while (t_start_point_idx < t_num_link_start_voxel) && num_unvisited_points > 0
    % Find the starting voxel list index in the voxel list - for links
    for t_start_point_idx = t_start_search_idx : t_num_link_start_voxel
        t_num_cc_voxel = 0;   
        t_current_id = l_link_start_voxel_idx(t_start_point_idx);
        if voxel_unvisited_Q(t_current_id)
            t_start_search_idx = t_start_point_idx + 1;
            keep_tracking = true;
            break;
        end
    end
    while keep_tracking
        keep_tracking = false;
        % Add the current voxel to the connected component voxel list 
        t_num_cc_voxel = t_num_cc_voxel + 1;
        % MATLAB can extend the length of the list automatically if
        % t_num_cc_voxel is larger than the initialized length. 
        t_link_idx(t_num_cc_voxel) = t_current_id;
        voxel_unvisited_Q(t_current_id) = false;
        % Get the neighbors of the current voxel and pick the ONE hasn't
        % been visited. 
        t_neighbor_idx = voxel_neighbor_idx_list(:, t_current_id);
        for tmp_idx = 1 : 2
            tmp_id = t_neighbor_idx(tmp_idx);
            if tmp_id > 0
                if voxel_unvisited_Q(tmp_id)
                    keep_tracking = true;
                    t_current_id = tmp_id;
                    break;
                end
            else
                break;
            end
        end
    end
    if t_num_cc_voxel > 0
        t_num_cc = t_num_cc + 1;
        num_unvisited_points = num_unvisited_points - t_num_cc_voxel;
        link_cc.PixelIdxList{t_num_cc} = voxel_list(t_link_idx(1:t_num_cc_voxel));
    end
end
link_cc.PixelIdxList = link_cc.PixelIdxList(1:t_num_cc)';
link_cc.NumObjects = t_num_cc;

if nargout > 1
    if nargout == 2
        varargout{1} = num_unvisited_points;
    elseif nargout == 3
        varargout{1} = num_unvisited_points;
        varargout{2} = voxel_unvisited_Q;
    end
end
end