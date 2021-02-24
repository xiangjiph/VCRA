function nf = fun_graph_get_node_features(cc_ind, vessel_image, vessel_mask_dt, image_size, compute_feature)
% fun_graph_get_node_features computes the node features including average
% intensity, average radius, center of mass, radial node density
% distribution, local neighbor distance and label, nearest neighbor
% distance and label. 
% Input: 
%   cc_ind: cell array, each cell caontains a numerical vector that is the
%   linear indices of the node voxel
%   vessel_image: 3D numerical array
%   vessel_mask_dt: 3D single array, distance transform of the vessel mask
%   image_size: size of the numerical array that cc_ind is in. Used for
%   converting cc_ind to the 3D location in the array
%   compute_feature: cell array of strings 'geometry', 'int', and 'dt'. 
% Ouput:
%   nf: structure with fields:
%       dt_mean: numerical vector, mean of the radius of the node voxels
%       int_mean: numerical vector, mean of the intenisty of the node voxels
%       nearest_node_dist: numerical vector, distnace to the nearest
%       neighbor node
%       nearest_node_label: numerical vector, label of the nearest neighbor
%       node
%       center_of_mass: N-by-3 numerical vector, the center of mass
%       location of the node in the array
%       radial_node_count, radial_node_count_edge: histogram count
%       of the number of node in the radial direction 
%       neighbor_node_label    
%       radial_node_density: radial_node_count normalized by the volume of
%       the radial sphere shell
%       neighbor_node_label, neighbor_node_dist: label and distnace of
%       nodes within a threshold distant from the node
%
% Implemented by Xiang Ji on 02/24/2019

%% Phase input
if nargin < 4
    image_size = [];
    compute_feature = {'geometry', 'int', 'dt'};
elseif nargin < 5
    compute_feature = {'geometry', 'int', 'dt'};
end
num_node = numel(cc_ind);

if isempty(image_size)
    if ~isempty(vessel_image)
        image_size = size(vessel_image);
    elseif ~isempty(vessel_mask_dt)
        image_size = size(vessel_mask_dt);
    else
        error('Both vessel_image and vessel_mask_dt are empty');
    end
end

if any(ismember(compute_feature, 'geometry'))
    compute_geometryQ = true;
else
    compute_geometryQ = false;
end
if any(ismember(compute_feature, 'int'))
    compute_intQ = true;
    if isempty(vessel_image)
        error('Vessel image is empty');
    end
else
    compute_intQ = false;
end
if any(ismember(compute_feature, 'dt'))
    compute_dtQ = true;
    if isempty(vessel_mask_dt)
        error('Vessel_mask_dt is empty');
    end
else
    compute_dtQ = false;
end
%% Computation
nf = struct;
[nf.dt_mean, nf.int_mean, nf.nearest_node_dist, nf.nearest_node_label] = deal(zeros(num_node, 1));
nf.center_of_mass = zeros(3, num_node);
for iter_node = 1 : num_node
    tmp_ind = cc_ind{iter_node};
    tmp_sub = fun_ind2sub(image_size, tmp_ind);
    if compute_geometryQ
        if isvector(tmp_sub)
            nf.center_of_mass(:, iter_node) = tmp_sub;
        else
            nf.center_of_mass(:, iter_node) = mean(tmp_sub, 1);
        end
    end
    
    if compute_intQ
        nf.int_mean(iter_node) = mean(vessel_image(tmp_ind));
    end
    if compute_dtQ
        nf.dt_mean(iter_node) = mean(vessel_mask_dt(tmp_ind));
    end
end
% The following features is not very useful for graph refinement. They are
% useful for network analysis.
% radial_binning_edge = [0:10:100, 125:25:500, 550:50:1000, 1100:100:5000];
% record_neighbor_label_r = 50;
% if compute_geometryQ
%     nf.center_of_mass = nf.center_of_mass';
%     
%     nf.distance_to_boundary = min(min(nf.center_of_mass, image_size - nf.center_of_mass + 1), [], 2);
%     [nf.radial_node_count, nf.radial_node_count_edge, nf.neighbor_node_label, ...
%         nf.neighbor_node_dist, nf.radial_node_density] = deal(cell(num_node, 1));
%     node_pdist = squareform(pdist(nf.center_of_mass));
%     for iter_node = 1 : num_node
%         tmp_node_neighbor_dist = node_pdist(:, iter_node);
%         tmp_valid_neighbor_label = find(tmp_node_neighbor_dist <= nf.distance_to_boundary(iter_node) & tmp_node_neighbor_dist > 0);
%         
%         tmp_valid_neighbor_dist = tmp_node_neighbor_dist(tmp_valid_neighbor_label);
%         tmp_bin = radial_binning_edge(radial_binning_edge < nf.distance_to_boundary(iter_node));
%         if numel(tmp_bin) > 1
%             [tmp_count, tmp_edge] = histcounts(tmp_valid_neighbor_dist, tmp_bin);
%             tmp_radial_layer_vol = 4 * pi / 3 .* (tmp_edge(2:end).^3 - tmp_edge(1:end-1).^3);
%             nf.radial_node_count{iter_node} = tmp_count;
%             nf.radial_node_count_edge{iter_node} = tmp_edge;
%             nf.radial_node_density{iter_node} = tmp_count ./ tmp_radial_layer_vol;
%         end
%         tmp_record_neighbor_Q = (tmp_valid_neighbor_dist <= record_neighbor_label_r);
%         if any(tmp_record_neighbor_Q)
%             tmp_record_neighbor_label = tmp_valid_neighbor_label(tmp_record_neighbor_Q);
%             tmp_record_neighbor_dt = tmp_valid_neighbor_dist(tmp_record_neighbor_Q);
%             % Sort 
%             [tmp_record_neighbor_dt, tmp_idx] = sort(tmp_record_neighbor_dt, 'ascend');
%             tmp_record_neighbor_label = tmp_record_neighbor_label(tmp_idx);
%             nf.nearest_node_dist(iter_node) = tmp_record_neighbor_dt(1);
%             nf.nearest_node_label(iter_node) = tmp_record_neighbor_label(1);
%             nf.neighbor_node_label{iter_node} = tmp_record_neighbor_label;
%             nf.neighbor_node_dist{iter_node} = tmp_record_neighbor_dt;
%         else
%             nf.nearest_node_dist(iter_node) = nan;
%             nf.nearest_node_label(iter_node) = nan;
%         end
%     end
% end
end