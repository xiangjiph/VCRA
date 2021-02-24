function wholebrain_stat_str = fun_analysis_load_whole_brain_features_v1(grid_info, recon_mask_name, layer_list)


%% Initialization
DataManager = FileManager;
dataset_name = grid_info.dataset_name;
stack = grid_info.stack;
if nargin < 3
    layer_list = grid_info.layer;
end
num_layer_to_process = numel(layer_list);
num_core = 8;
%%
wb_grid_size = [grid_info.grid_size(1:2), num_layer_to_process];
[wb_link_feature_cell, wb_node_feature_cell, wb_mask_1_proj_cell, wb_mask_2_proj_cell, ...
    wb_mask_3_proj_cell] = deal(cell(wb_grid_size));
[wb_volume_ratio_cell, wb_mask_surf_area_2_vol_cell_mm2mm3] = deal(cell(num_layer_to_process, 1));
%%
pool_obj = gcp('nocreate');
delete(pool_obj);
parpool(num_core);
parfor iter_layer = 1 : num_layer_to_process
    tmp_tic = tic;
    % Initialization within layer
    grid_layer = layer_list(iter_layer);
    fprintf('Load layer %d\n', grid_layer);
    block_grid_sub = grid_info.bbox_grid_sub{grid_layer};
    num_valid_graph = grid_info.num_bbox_xy(grid_layer);
    [layer_link_feature_cell, layer_node_feature_cell, layer_mask_1_proj_cell, ...
        layer_mask_2_proj_cell, layer_mask_3_proj_cell] = deal(cell(grid_info.grid_size(1:2)));
    [layer_volume_ratio, layer_surf_area2vol_mm2pmm3] = deal(nan(grid_info.grid_size(1:2)));
    % Load graph in single plane into cell arrays
    for iter_block = 1 : num_valid_graph
        tmp_idx_1 = block_grid_sub(iter_block, 1);
        tmp_idx_2 = block_grid_sub(iter_block, 2);
        tmp_layer = block_grid_sub(iter_block, 3);
        try
            tmp_mask_str = DataManager.load_block_mask(dataset_name, stack, recon_mask_name, ...
                tmp_idx_1, tmp_idx_2, tmp_layer);
            tmp_mask_str.block_volume_mm3 = prod(tmp_mask_str.block_size) / 1e9;
            layer_mask_1_proj_cell{tmp_idx_1, tmp_idx_2} = tmp_mask_str.max_proj_1;
            layer_mask_2_proj_cell{tmp_idx_1, tmp_idx_2} = tmp_mask_str.max_proj_2;
            layer_mask_3_proj_cell{tmp_idx_1, tmp_idx_2} = tmp_mask_str.max_proj_3;
            layer_link_feature_cell{tmp_idx_1, tmp_idx_2} = tmp_mask_str.link.features;
            layer_node_feature_cell{tmp_idx_1, tmp_idx_2} = tmp_mask_str.node.features;
            
            layer_volume_ratio(tmp_idx_1, tmp_idx_2) = tmp_mask_str.mask_volume_ratio;
            % Normalized to mm^2 / mm^3
            layer_surf_area2vol_mm2pmm3(tmp_idx_1, tmp_idx_2) = tmp_mask_str.mask_surface_area / (1e6 * tmp_mask_str.block_volume_mm3) ;
        catch ME
            fprintf('Unable to load reconstructed mask (%d, %d, %d)\n', tmp_idx_1, tmp_idx_2, tmp_layer);
            fprintf('Error message: %s\n', ME.message);
            continue;
        end
    end
    wb_link_feature_cell(:, :, iter_layer) = layer_link_feature_cell;
    wb_node_feature_cell(:, :, iter_layer) = layer_node_feature_cell;
    wb_mask_1_proj_cell(:, :, iter_layer) = layer_mask_1_proj_cell;
    wb_mask_2_proj_cell(:, :, iter_layer) = layer_mask_2_proj_cell;
    wb_mask_3_proj_cell(:, :, iter_layer) = layer_mask_3_proj_cell;
    wb_volume_ratio_cell{iter_layer} = layer_volume_ratio;
    wb_mask_surf_area_2_vol_cell_mm2mm3{iter_layer} = layer_surf_area2vol_mm2pmm3;
    fprintf('Finish processing layer %d. Elapsed time is %f seconds\n', grid_layer, toc(tmp_tic));
end
wb_volume_ratio = cat(3, wb_volume_ratio_cell{:});
wb_total_surf_area = cat(3, wb_mask_surf_area_2_vol_cell_mm2mm3{:});
%% Save to structure
wholebrain_stat_str = struct;
wholebrain_stat_str.grid_size = wb_grid_size;
wholebrain_stat_str.link_features = wb_link_feature_cell;
wholebrain_stat_str.node_features = wb_node_feature_cell;
wholebrain_stat_str.mask_max_proj_1 = wb_mask_1_proj_cell;
wholebrain_stat_str.mask_max_proj_2 = wb_mask_2_proj_cell;
wholebrain_stat_str.mask_max_proj_3 = wb_mask_3_proj_cell;
wholebrain_stat_str.volume_ratio = wb_volume_ratio;
wholebrain_stat_str.surface_area_mm2mm3 = wb_total_surf_area;
end