function dtf = fun_analysis_reonstruction_space_properties(mask, mask_labeled, save_cc_voxel_infoQ)
% fun_analysis_reonstruction_space_properties computes the distance between
% tissue to its nearest vessels, the position of the voxels that are
% closest to each node and link in the reconstruction mask and perform
% basic analysis. 
% Input:
%   mask: 3D logical array, reconstructed vessel mask 
%   mask_labeled: 3D int array, reconstructed vessel mask labeled with the
%   node and link label from whom the mask is reconstructed from. Positive
%   integer voxel value is the link label the voxel belong to. Negative
%   integer voxel value is the opposite of the node label that the voxel
%   belong to. 
%   save_cc_voxel_info: logical scalar. If true, returned structure will
%   include the nearest tissue voxel position cell array and the their
%   distanctance transform value
% Output:
%   dtf: distance transform features. Structure with fields 
%
% Implemented by Xiang Ji on 02/25/2019

if nargin < 3
    save_cc_voxel_infoQ = false;
end
dtf = struct;
dtf.Mask_volume= nnz(mask);
dtf.Block_volume = numel(mask);
image_size = size(mask);

% Compute distance between the tissue to the nearest vessel
[tissue_space_dt, tissue_nearest_mask_ind] = bwdist(mask);
valid_tissue_voxel_ind = find(tissue_space_dt > 0);
if dtf.Block_volume  < intmax('uint16')
    valid_tissue_voxel_ind = uint16(valid_tissue_voxel_ind);
elseif dtf.Block_volume  < intmax('uint32')
    valid_tissue_voxel_ind = uint32(valid_tissue_voxel_ind);
end
% The mask voxel indices (before cropping) that are closest to the space point
% valid_nearest_mask_ind_list = tissue_nearest_mask_ind(valid_tissue_voxel_ind);
valid_tissue_space_dt_list = tissue_space_dt(valid_tissue_voxel_ind);
clearvars tissue_space_dt mask

% The nearest label is from the graph that reconstruct the labeled matrix,
% not the graph that use for internal analysis
% [space_idx_list, nearest_cc_label] = fun_bin_data_to_idx_list(mask_labeled(valid_nearest_mask_ind_list));
[space_idx_list, nearest_cc_label] = fun_bin_data_to_idx_list(mask_labeled(tissue_nearest_mask_ind(valid_tissue_voxel_ind)));
clearvars tissue_nearest_mask_ind mask_labeled
% Convert the list index to voxel indices in the block 
num_label = numel(nearest_cc_label);
[space_ind_list, space_dt_list] = deal(cell(num_label, 1));
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
%     completely_valid_Q(iter_label) = ~tmp_on_boundaryQ;
    if ~tmp_on_boundaryQ
        tmp_dt = valid_tissue_space_dt_list(tmp_list_idx);
%     % Distance to the edge of the mask 
%     tmp_dt_to_edge = min(min(tmp_sub, image_size - tmp_sub), [], 2);
%     tmp_validQ = tmp_dt_to_edge > tmp_dt;
%     completely_valid_Q(iter_label) = all(tmp_validQ);
        if save_cc_voxel_infoQ 
            space_ind_list{iter_label} = tmp_ind;
        end
        space_dt_list{iter_label} = tmp_dt;
        space_vol(iter_label) = numel(tmp_dt);
        space_dt_mean(iter_label) = mean(tmp_dt);
        space_dt_max(iter_label) = max(tmp_dt);
        space_dt_median(iter_label) = median(tmp_dt);
    end
end
% clearvars valid_nearest_mask_ind_list valid_tissue_space_dt_list
dtf.global_stat = fun_analysis_get_basic_statistics(cat(1, space_dt_list{:}));

closest_to_link_Q = (nearest_cc_label > 0);
closest_to_node_Q = (nearest_cc_label < 0);

dtf.closest_link_label = nearest_cc_label(closest_to_link_Q);
dtf.closest_node_label = - nearest_cc_label(closest_to_node_Q);
dtf.node_nearest_tissue_dt_mean = space_dt_mean(closest_to_node_Q);
dtf.link_nearest_tissue_dt_mean = space_dt_mean(closest_to_link_Q);

dtf.node_nearest_tissue_dt_max = space_dt_max(closest_to_node_Q);
dtf.link_nearest_tissue_dt_max = space_dt_max(closest_to_link_Q);

dtf.node_nearest_tissue_dt_median = space_dt_median(closest_to_node_Q);
dtf.link_nearest_tissue_dt_median = space_dt_median(closest_to_link_Q);

dtf.node_nearest_tissue_volume = space_vol(closest_to_node_Q);
dtf.link_nearest_tissue_volume = space_vol(closest_to_link_Q);

if save_cc_voxel_infoQ
    dtf.node_nearest_tissue_ind = space_ind_list(closest_to_node_Q);
    dtf.link_nearest_tissue_ind = space_ind_list(closest_to_link_Q);
    
    dtf.node_nearest_tissue_dt = space_dt_list(closest_to_node_Q);
    dtf.link_nearest_tissue_dt = space_dt_list(closest_to_link_Q);
end
end