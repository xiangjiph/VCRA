function wholebrain_stat_str = fun_analysis_load_whole_brain_features(grid_info, recon_mask_name, layer_list)


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
% Anisotropy
[wb_anisotropy_all_vw_fa, wb_anisotropy_cap_vw_fa, ...
    wb_anisotropy_all_vw_fa_z, wb_anisotropy_cap_vw_fa_z, ...
    wb_anisotropy_all_vw_min2max_z, wb_anisotropy_cap_vw_min2max_z, ...
    wb_anisotropy_all_vw_svd1, wb_anisotropy_cap_vw_svd1, ...
    wb_anisotropy_all_vw_svd1_z, wb_anisotropy_cap_vw_svd1_z] = deal(cell(num_layer_to_process, 1));
[wb_anisotropy_all_vw_vec, wb_anisotropy_cap_vw_vec] = deal(cell(num_layer_to_process, 1));
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
    % Record volume-weighted anisotropy information
    [layer_anisotropy_all_vw_fa, layer_anisotropy_cap_vw_fa, ...
        layer_anisotropy_all_vw_fa_z, layer_anisotropy_cap_vw_fa_z,...
        layer_anisotropy_all_vw_min2max_z, layer_anisotropy_cap_vw_min2max_z, ...
        layer_anisotropy_all_vw_svd1, layer_anisotropy_cap_vw_svd1, ...
        layer_anisotropy_all_vw_svd1_z, layer_anisotropy_cap_vw_svd1_z] = deal(nan(grid_info.grid_size(1:2)));
    [layer_anisotropy_all_vw_vec, layer_anisotropy_cap_vw_vec] = deal(nan([3, grid_info.grid_size(1:2)]));
    
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
            % Record anisotropy
            if isfield(tmp_mask_str, 'anisotropy_all_vw') && isfield(tmp_mask_str.anisotropy_all_vw, 'min2max_z') ...
                    && ~isempty(tmp_mask_str.anisotropy_all_vw.min2max_z)
                layer_anisotropy_all_vw_fa(tmp_idx_1, tmp_idx_2) = tmp_mask_str.anisotropy_all_vw.fractional_anisotropy;
                layer_anisotropy_all_vw_fa_z(tmp_idx_1, tmp_idx_2) = tmp_mask_str.anisotropy_all_vw.fa_z;
                layer_anisotropy_all_vw_min2max_z(tmp_idx_1, tmp_idx_2) = tmp_mask_str.anisotropy_all_vw.min2max_z;
                layer_anisotropy_all_vw_svd1_z(tmp_idx_1, tmp_idx_2) = tmp_mask_str.anisotropy_all_vw.svd_1_z;
                layer_anisotropy_all_vw_svd1(tmp_idx_1, tmp_idx_2) = tmp_mask_str.anisotropy_all_vw.svd_value_ratio(1);
                layer_anisotropy_all_vw_vec(:, tmp_idx_1, tmp_idx_2) = tmp_mask_str.anisotropy_all_vw.svd_max_vec;
            end            
            if isfield(tmp_mask_str, 'anisotropy_capillary_vw') && isfield(tmp_mask_str.anisotropy_capillary_vw, 'min2max_z') ...
                    && ~isempty(tmp_mask_str.anisotropy_capillary_vw.min2max_z)
                layer_anisotropy_cap_vw_fa(tmp_idx_1, tmp_idx_2) = tmp_mask_str.anisotropy_capillary_vw.fractional_anisotropy;
                layer_anisotropy_cap_vw_fa_z(tmp_idx_1, tmp_idx_2) = tmp_mask_str.anisotropy_capillary_vw.fa_z;
                layer_anisotropy_cap_vw_min2max_z(tmp_idx_1, tmp_idx_2) = tmp_mask_str.anisotropy_capillary_vw.min2max_z;
                layer_anisotropy_cap_vw_svd1(tmp_idx_1, tmp_idx_2) = tmp_mask_str.anisotropy_capillary_vw.svd_value_ratio(1);
                layer_anisotropy_cap_vw_svd1_z(tmp_idx_1, tmp_idx_2) = tmp_mask_str.anisotropy_capillary_vw.svd_1_z;
                layer_anisotropy_cap_vw_vec(:, tmp_idx_1, tmp_idx_2) = tmp_mask_str.anisotropy_capillary_vw.svd_max_vec;                
            end            
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
    % Anisotropy
    wb_anisotropy_all_vw_fa{iter_layer} = layer_anisotropy_all_vw_fa;
    wb_anisotropy_all_vw_fa_z{iter_layer} = layer_anisotropy_all_vw_fa_z;
    wb_anisotropy_all_vw_min2max_z{iter_layer} = layer_anisotropy_all_vw_min2max_z;
    wb_anisotropy_all_vw_svd1{iter_layer} = layer_anisotropy_all_vw_svd1;
    wb_anisotropy_all_vw_svd1_z{iter_layer} = layer_anisotropy_all_vw_svd1_z;
    wb_anisotropy_all_vw_vec{iter_layer} = layer_anisotropy_all_vw_vec;
    
    wb_anisotropy_cap_vw_fa{iter_layer} = layer_anisotropy_cap_vw_fa;
    wb_anisotropy_cap_vw_fa_z{iter_layer} = layer_anisotropy_cap_vw_fa_z;
    wb_anisotropy_cap_vw_min2max_z{iter_layer} = layer_anisotropy_cap_vw_min2max_z;
    wb_anisotropy_cap_vw_svd1{iter_layer} = layer_anisotropy_cap_vw_svd1;
    wb_anisotropy_cap_vw_svd1_z{iter_layer} = layer_anisotropy_cap_vw_svd1_z;
    wb_anisotropy_cap_vw_vec{iter_layer} = layer_anisotropy_cap_vw_vec;
    fprintf('Finish processing layer %d. Elapsed time is %f seconds\n', grid_layer, toc(tmp_tic));
end
%% Save to structure
wholebrain_stat_str = struct;
wholebrain_stat_str.grid_size = wb_grid_size;
wholebrain_stat_str.link_features = wb_link_feature_cell;
wholebrain_stat_str.node_features = wb_node_feature_cell;
wholebrain_stat_str.mask_max_proj_1 = wb_mask_1_proj_cell;
wholebrain_stat_str.mask_max_proj_2 = wb_mask_2_proj_cell;
wholebrain_stat_str.mask_max_proj_3 = wb_mask_3_proj_cell;
wholebrain_stat_str.volume_ratio = cat(3, wb_volume_ratio_cell{:});
wholebrain_stat_str.surface_area_mm2mm3 = cat(3, wb_mask_surf_area_2_vol_cell_mm2mm3{:});
% Anisotropy - volume weighted
wholebrain_stat_str.anisotropy_all_vw_fa = cat(3, wb_anisotropy_all_vw_fa{:});
wholebrain_stat_str.anisotropy_all_vw_fa_z = cat(3, wb_anisotropy_all_vw_fa_z{:});
wholebrain_stat_str.anisotropy_all_vw_min2max_z = cat(3, wb_anisotropy_all_vw_min2max_z{:});
wholebrain_stat_str.anisotropy_all_vw_svd1 = cat(3, wb_anisotropy_all_vw_svd1{:});
wholebrain_stat_str.anisotropy_all_vw_svd1_z = cat(3, wb_anisotropy_all_vw_svd1_z{:});

wholebrain_stat_str.anisotropy_cap_vw_fa = cat(3, wb_anisotropy_cap_vw_fa{:});
wholebrain_stat_str.anisotropy_cap_vw_fa_z = cat(3, wb_anisotropy_cap_vw_fa_z{:});
wholebrain_stat_str.anisotropy_cap_vw_min2max_z = cat(3, wb_anisotropy_cap_vw_min2max_z{:});
wholebrain_stat_str.anisotropy_cap_vw_svd1 = cat(3, wb_anisotropy_cap_vw_svd1{:});
wholebrain_stat_str.anisotropy_cap_vw_svd1_z = cat(3, wb_anisotropy_cap_vw_svd1_z{:});

wholebrain_stat_str.anisotropy_all_vw_vec = cat(4, wb_anisotropy_all_vw_vec{:});
wholebrain_stat_str.anisotropy_cap_vw_vec = cat(4, wb_anisotropy_cap_vw_vec{:});
end