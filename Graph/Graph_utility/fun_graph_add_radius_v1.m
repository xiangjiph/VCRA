function graph = fun_graph_add_radius(graph, distance_transform, add_type, min_1_Q)
% fun_graph_add_radius adds radius from the distance transform to the graph
% structure. 
% Input:
%   graph: structure output by fun_skeleton_to_graph
%   distance_transform: distance transform of the mask that the graph was
%   derived from. 
%   add_type: determine the type of the added field
% Output: 
%   graph: graph structure with radius of each node, link, endpoint and
%   isolated points added. 
% Note: 
%   For later graph stitching purpose, the radia are stored as cell arrays
%   for each connected components of nodes and links. 
if nargin < 3
    add_type = 'sparse';
    min_1_Q = true;
elseif nargin < 4
    min_1_Q = true;
end

switch add_type
    case 'cell'
        % Add radius for node
        graph.node.radius = cell(graph.node.num_cc, 1);
        for tmp_idx = 1 : graph.node.num_cc
            if min_1_Q
                graph.node.radius{tmp_idx} = min(1, distance_transform(graph.node.cc_ind{tmp_idx}));
            else
                graph.node.radius{tmp_idx} = distance_transform(graph.node.cc_ind{tmp_idx});
            end
        end
        % Add radius for link
        graph.link.radius = cell(graph.link.num_cc, 1);
        for tmp_idx = 1 : graph.link.num_cc
            if min_1_Q
                graph.link.radius{tmp_idx} = min(1, distance_transform(graph.link.cc_ind{tmp_idx}));
            else
                graph.link.radius{tmp_idx} = distance_transform(graph.link.cc_ind{tmp_idx});
            end
        end

        % Add radius for isolated points
        if isfield(graph, 'isopoint')
            graph.isopoint.radius = zeros(graph.isopoint.num_voxel, 1);
            for tmp_idx = 1 : graph.isopoint.num_voxel
                if min_1_Q
                    graph.isopoint.radius(tmp_idx) = min(1, distance_transform(graph.isopoint.pos_ind(tmp_idx)));
                else
                    graph.isopoint.radius(tmp_idx) = distance_transform(graph.isopoint.pos_ind(tmp_idx));
                end
            end
        end

        % Add radius for end points
        if isfield(graph, 'endpoint')
            graph.endpoint.radius = zeros(graph.endpoint.num_voxel, 1);
            for tmp_idx = 1 : graph.endpoint.num_voxel
                if min_1_Q
                    graph.endpoint.radius(tmp_idx) = min(1, distance_transform(graph.endpoint.pos_ind(tmp_idx)));
                else
                    graph.endpoint.radius(tmp_idx) = distance_transform(graph.endpoint.pos_ind(tmp_idx));
                end
            end
        end
    case 'vector'
        if min_1_Q
            graph.node.radius = min(1, distance_transform(graph.node.pos_ind));
            graph.link.radius = min(1, distance_transform(graph.link.pos_ind));
            graph.isopoint.radius = min(1, distance_transform(graph.isopoint.pos_ind));
            graph.endpoint.radius = min(1, distance_transform(graph.endpoint.pos_ind));
        else
            graph.node.radius = distance_transform(graph.node.pos_ind);
            graph.link.radius = distance_transform(graph.link.pos_ind);
            graph.isopoint.radius = distance_transform(graph.isopoint.pos_ind);
            graph.endpoint.radius = distance_transform(graph.endpoint.pos_ind);
        end
    case 'sparse'
        % Get the list of all the voxels in the skeleton: 
        % The link voxel contains the end point voxels. 
        % Assuming *.pos_ind = cat(1, *.cc_ind{:})
        voxel_list = cat(1, graph.node.pos_ind, graph.link.pos_ind, ...
            graph.isopoint.pos_ind(:));
        if min_1_Q
            graph.radius = sparse(voxel_list, ones(numel(voxel_list),1), ...
                min(1, double(distance_transform(voxel_list))), numel(distance_transform),1);
        else
            graph.radius = sparse(voxel_list, ones(numel(voxel_list),1), ...
                double(distance_transform(voxel_list)), numel(distance_transform),1);
        end
        
end
end