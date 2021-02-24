function graph = fun_graph_add_radius(graph, mask_dt, min_radius)
% fun_graph_add_radius adds radius from the distance transform to the graph
% structure. 
% Input:
%   graph: structure output by fun_skeleton_to_graph
%   mask_dt: distance transform of the mask that the graph was
%   derived from. 
%   min_radius: If greater than 0, set the minimum radius of the voxel to be this value 
% Output: 
%   graph: graph structure with radius of each node, link
if nargin < 3
    min_radius = 0.25;
end
if islogical(min_radius)
    min_radius = double(min_radius);
end

voxel_ind = cat(1, graph.node.cc_ind, graph.link.cc_ind);
if ~isempty(graph.isopoint.pos_ind)
    voxel_ind = cat(1, voxel_ind, num2cell(graph.isopoint.pos_ind));
end
if isfield(graph, 'isoloop') && ~isempty(graph.isoloop.cc_ind)
    voxel_ind = cat(1, voxel_ind, graph.isoloop.cc_ind);
end

num_cc = numel(voxel_ind);
voxel_rad = cell(num_cc, 1);
if issparse(mask_dt)
    dt_is_sparse_Q = true;
else
    dt_is_sparse_Q = false;
end

for iter_cc = 1 : num_cc
    tmp_ind = voxel_ind{iter_cc};
    if dt_is_sparse_Q
        tmp_dt = full(mask_dt(tmp_ind));
    else
        tmp_dt = mask_dt(tmp_ind);
    end
    tmp_dt_nonzero_Q = tmp_dt>0;
    if any(tmp_dt_nonzero_Q)
        tmp_median = median(tmp_dt(tmp_dt_nonzero_Q));
        % If the radius of the voxel is 0, replace it by the median radius
        % of this segment
        tmp_dt(~tmp_dt_nonzero_Q) = tmp_median;
        voxel_rad{iter_cc} = tmp_dt;
    else
        tmp_num_voxel = numel(tmp_ind);
        fprintf('Exist cc with %d voxels that is completely outside the distance transform mask\n', tmp_num_voxel);
        voxel_rad{iter_cc} = repelem(0, tmp_num_voxel, 1);
    end
end
voxel_rad = cat(1, voxel_rad{:});
voxel_ind = cat(1, voxel_ind{:});
if min_radius
%     if min_radius > 0, set the minimum voxel radius to 0
    graph.radius = sparse(voxel_ind, ones(numel(voxel_ind),1), ...
        max(min_radius, double(voxel_rad)), numel(mask_dt),1);
else
    graph.radius = sparse(voxel_ind, ones(numel(voxel_ind),1), ...
        double(voxel_rad), numel(mask_dt),1);
end
        
end