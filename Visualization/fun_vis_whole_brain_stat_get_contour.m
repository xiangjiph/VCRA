function vis_contour_str = fun_vis_whole_brain_stat_get_contour(vis_layer_bbox_mmxx, registration_vis_str,...
    vis_projection_dim, vis_plan_dim, vis_downsample_rate)

if isempty(vis_layer_bbox_mmxx)
    vis_contour_str = [];
else    
    vis_contour_plane_med = vis_layer_bbox_mmxx(1, [vis_projection_dim, vis_projection_dim + 3]);
    vis_contour_plane_med = round((vis_contour_plane_med(1) + vis_contour_plane_med(2))/2 ...
        ./ registration_vis_str.voxel_size(vis_projection_dim));
    vis_contour_plane_labeled_array = registration_vis_str.combined_regional_mask(:, :, ...
        vis_contour_plane_med);
    vis_hemisphere_plane_array = registration_vis_str.registered_hemisphere_mask(:, :, ...
        vis_contour_plane_med);
    
    vis_unique_label = unique(vis_contour_plane_labeled_array, 'sorted');
    vis_unique_label = vis_unique_label(vis_unique_label > 0);
    num_vis_lable = numel(vis_unique_label);
    vis_contour_str = struct;
    vis_contour_str.sub_cell = cell(2, num_vis_lable);
    for iter_label = 1 : num_vis_lable
        tmp_contour_plane_labeled_array = (vis_contour_plane_labeled_array == vis_unique_label(iter_label));
        tmp_right_hemisphere_contour_labeled_array = (tmp_contour_plane_labeled_array & vis_hemisphere_plane_array);
        tmp_left_hemisphere_contour_labeled_array = (tmp_contour_plane_labeled_array & ~vis_hemisphere_plane_array);
        
        vis_contour_str.sub_cell{1, iter_label} = cellfun(@(x) (registration_vis_str.voxel_size(vis_plan_dim) ...
            ./ vis_downsample_rate) .* x, bwboundaries(tmp_right_hemisphere_contour_labeled_array, 8), ...
            'UniformOutput', false);
        vis_contour_str.sub_cell{2, iter_label} = cellfun(@(x) (registration_vis_str.voxel_size(vis_plan_dim) ...
            ./ vis_downsample_rate) .* x, bwboundaries(tmp_left_hemisphere_contour_labeled_array, 8), ...
            'UniformOutput', false);
    end
    is_not_empty_cell_Q = ~cellfun(@isempty, vis_contour_str.sub_cell);
    vis_num_side = sum(is_not_empty_cell_Q, 1);
    vis_unique_label = repelem(vis_unique_label, vis_num_side, 1);
    vis_contour_str.sub_cell = vis_contour_str.sub_cell(is_not_empty_cell_Q);
    
    vis_contour_str.structure_name = cellfun(@(x) registration_vis_str.label_2_name(x), ...
        num2cell(vis_unique_label), 'UniformOutput', false);
    vis_contour_str.structure_abbrv = cellfun(@(x) registration_vis_str.label_2_acronym(x), ...
        num2cell(vis_unique_label), 'UniformOutput', false);
end
end