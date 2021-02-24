function noncapillary_DT_str = fun_analysis_reconstruction_noncapillary_space_properties(input_graph, capillary_max_radius, scale, recon_max_error_rate)
% fun_analysis_reconstruction_noncapillary_space_properties does the
% following: 
% 1. Select the vessel segments in the graph that is of median radius
% greater than capillary_max_radius, reconstruct the non-capillary mask and
% calculate the distance transform field outside the non-capillary mask. 
% 2. Compute features of the extra-noncapillary space 
% 3. Compute the number of capillary in the nearest tissue volume
% associated with each reconstructed noncapillary segment. 
%% Default value
if nargin < 3
    scale = 0.5;
    recon_max_error_rate = 0.1;
elseif nargin < 4
    recon_max_error_rate = 0.1;
end
%% Initialization
noncapillary_DT_str = struct;
[noncapillary_DT_str.noncapillary_nearest_tissue_dt_mean, ...
    noncapillary_DT_str.noncapillary_nearest_tissue_dt_max, ...
    noncapillary_DT_str.noncapillary_nearest_tissue_dt_median, ...
    noncapillary_DT_str.noncapillary_nearest_tissue_dt_volume, ...
    noncapillary_DT_str.noncapillary_nearest_capillary_num_vxl, ...
    noncapillary_DT_str.dist_ep1_to_nearest_noncapillary, ...
    noncapillary_DT_str.dist_ep2_to_nearest_noncapillary, ...
    noncapillary_DT_str.dist_to_nearest_noncapillary_mean, ...
    noncapillary_DT_str.dist_to_nearest_noncapillary_median] = deal(nan(input_graph.link.num_cc, 1));
if input_graph.link.num_cc == 1
    noncapillary_DT_str = struct2table(noncapillary_DT_str, 'AsArray', true);
else
    noncapillary_DT_str = struct2table(noncapillary_DT_str);
end
if input_graph.link.num_cc == 0
    return;
end
%% Select large vessels for reconstruction 
all_vessel_radius = fun_analysis_get_cc_median_radius(input_graph.link.cc_ind, input_graph.radius);
is_large_vessel_Q = all_vessel_radius > capillary_max_radius;
if ~any(is_large_vessel_Q)
    return;
end
large_vessel_link_label = find(is_large_vessel_Q);

large_vessel_link_cc = input_graph.link.cc_ind(large_vessel_link_label);

large_vessel_node_label = input_graph.link.connected_node_label(large_vessel_link_label, :);
large_vessel_node_label = unique(large_vessel_node_label(large_vessel_node_label > 0));
large_vessel_node_cc = input_graph.node.cc_ind(large_vessel_node_label);

large_vessel_recon_cc = cat(1, large_vessel_node_cc, large_vessel_link_cc);

large_vessel_recon_label_node = repelem(large_vessel_node_label, cellfun(@numel, large_vessel_node_cc));
large_vessel_recon_label_link = repelem(large_vessel_link_label, cellfun(@numel, large_vessel_link_cc));
if isrow(large_vessel_recon_label_node)
    large_vessel_recon_label_node = large_vessel_recon_label_node.';
end
if isrow(large_vessel_recon_label_link)
    large_vessel_recon_label_link = large_vessel_recon_label_link.';
end
large_vessel_recon_label = cat(1, - large_vessel_recon_label_node, large_vessel_recon_label_link);

large_vessel_recon_r = fun_graph_get_cc_radius(large_vessel_recon_cc, input_graph.radius, 'all');
large_vessel_recon_ind = cat(1, large_vessel_recon_cc{:});

mask_labeled = fun_skeleton_reconstruction_label_aprox(large_vessel_recon_ind, ...
    large_vessel_recon_r, large_vessel_recon_label, input_graph.num.mask_size, recon_max_error_rate);
%% Compute distance transform 
image_size = size(mask_labeled);
im_dims = nnz(image_size ~= 1);

if scale ~= 1
    mask_labeled = imresize3(mask_labeled, scale, 'Method', 'nearest');
end
% Compute distance between the tissue to the nearest vessel
[tissue_space_dt, tissue_nearest_mask_ind] = bwdist(mask_labeled~=0);
dt_size = size(tissue_space_dt);
%% Convert the capillary cc_ind to the resized coordinate
capillary_cc_ind = input_graph.link.cc_ind(~is_large_vessel_Q);
if scale ~= 1
    assert(all(abs(dt_size - image_size .* scale) < 1), 'Inconsistent array size before and after imresize3.');
    capillary_cc_ind = fun_analysis_get_scaled_cc_ind(capillary_cc_ind, ...
        image_size, scale);
end
% Get the distance properties of the capillary to the nearest large vessels
cap_to_noncap_dist_str = fun_analysis_get_capillary_to_large_vessel_DT_features(capillary_cc_ind, ...
    mask_labeled, tissue_space_dt .* (1/scale), tissue_nearest_mask_ind);

capillary_link_label = find(~is_large_vessel_Q);
noncapillary_DT_str.dist_ep1_to_nearest_noncapillary(capillary_link_label) = cap_to_noncap_dist_str.dist_to_nearest_noncapillary_ep1;
noncapillary_DT_str.dist_ep2_to_nearest_noncapillary(capillary_link_label) = cap_to_noncap_dist_str.dist_to_nearest_noncapillary_ep2;
noncapillary_DT_str.dist_to_nearest_noncapillary_mean(capillary_link_label) = cap_to_noncap_dist_str.dist_to_nearest_noncapillary_mean;
noncapillary_DT_str.dist_to_nearest_noncapillary_median(capillary_link_label) = cap_to_noncap_dist_str.dist_to_nearest_noncapillary_median;
% Get the length of capillary in the nearest tissue space for each
% noncapillary segment
[noncapillary_nearest_cap_num_vxl, noncapillary_label] = fun_bin_data_to_idx_list(cat(1, cap_to_noncap_dist_str.nearest_noncapillary_voxel_label{:}));
noncapillary_nearest_cap_num_vxl = cellfun(@numel, noncapillary_nearest_cap_num_vxl);
is_link_Q = noncapillary_label > 0;
noncapillary_label = noncapillary_label(is_link_Q);

noncapillary_DT_str.noncapillary_nearest_capillary_num_vxl(noncapillary_label) = noncapillary_nearest_cap_num_vxl(is_link_Q);
%%
num_voxel_in_mask = prod(image_size);
valid_tissue_voxel_ind = find(tissue_space_dt > 0);
if num_voxel_in_mask < intmax('uint16')
    valid_tissue_voxel_ind = uint16(valid_tissue_voxel_ind);
elseif num_voxel_in_mask  < intmax('uint32')
    valid_tissue_voxel_ind = uint32(valid_tissue_voxel_ind);
end
% The mask voxel indices (before cropping) that are closest to the space point
% valid_nearest_mask_ind_list = tissue_nearest_mask_ind(valid_tissue_voxel_ind);
valid_tissue_space_dt_list = tissue_space_dt(valid_tissue_voxel_ind);
clearvars tissue_space_dt;
% The nearest label is from the graph that reconstruct the labeled matrix,
% not the graph that use for internal analysis
[space_idx_list, nearest_cc_label] = fun_bin_data_to_idx_list(mask_labeled(tissue_nearest_mask_ind(valid_tissue_voxel_ind)));
clearvars tissue_nearest_mask_ind mask_labeled
% Convert the list index to voxel indices in the block 
% Only analysis the result of links. 
closest_to_link_Q = (nearest_cc_label > 0);
space_idx_list = space_idx_list(closest_to_link_Q);
nearest_cc_label = nearest_cc_label(closest_to_link_Q);

num_label = numel(nearest_cc_label);
[space_dt_list] = deal(cell(num_label, 1));
[space_dt_mean, space_dt_max, space_dt_median, space_vol] = deal(nan(num_label, 1));
% completely_valid_Q = false(num_label, 1);
for iter_label = 1 : num_label
    tmp_list_idx = space_idx_list{iter_label};
    tmp_ind = valid_tissue_voxel_ind(tmp_list_idx);
    [tmp_sub1, tmp_sub2, tmp_sub3] = ind2sub(image_size, tmp_ind);
    % Has voxel on the boundary:
    tmp_on_boundaryQ = any(tmp_sub1 == 1) || any(tmp_sub2 == 1) || any(tmp_sub3 == 1) ...
        || any(tmp_sub1 == image_size(1)) || any(tmp_sub2 == image_size(2)) || ...
        any(tmp_sub3 == image_size(3));
    if ~tmp_on_boundaryQ
        tmp_dt = valid_tissue_space_dt_list(tmp_list_idx);
%     % Distance to the edge of the mask 
%     tmp_dt_to_edge = min(min(tmp_sub, image_size - tmp_sub), [], 2);
%     tmp_validQ = tmp_dt_to_edge > tmp_dt;
%     completely_valid_Q(iter_label) = all(tmp_validQ);
        space_dt_list{iter_label} = tmp_dt;
        space_vol(iter_label) = numel(tmp_dt) ./ (scale ^ im_dims);
        space_dt_mean(iter_label) = mean(tmp_dt);
        space_dt_max(iter_label) = max(tmp_dt);
        space_dt_median(iter_label) = median(tmp_dt);
    end
end
clearvars valid_tissue_space_dt_list valid_tissue_voxel_ind
% Save
noncapillary_DT_str.noncapillary_nearest_tissue_dt_mean(nearest_cc_label) = space_dt_mean;
noncapillary_DT_str.noncapillary_nearest_tissue_dt_max(nearest_cc_label) = space_dt_max;
noncapillary_DT_str.noncapillary_nearest_tissue_dt_median(nearest_cc_label) = space_dt_median;
noncapillary_DT_str.noncapillary_nearest_tissue_dt_volume(nearest_cc_label) = space_vol;
end