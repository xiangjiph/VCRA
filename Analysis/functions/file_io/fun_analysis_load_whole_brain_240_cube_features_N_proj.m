function wholebrain_stat_str = fun_analysis_load_whole_brain_240_cube_features_N_proj(grid_info, recon_mask_name, layer_list)

%% Initialization
persistent DataManager
if isempty(DataManager)
    DataManager = FileManager;
end
dataset_name = grid_info.dataset_name;
stack = grid_info.stack;
if nargin < 3
    layer_list = grid_info.layer;
end
num_core = 8;
%% Simplify the parfor loading
load_bbox_grid_ind = cat(1, grid_info.bbox_grid_ind{layer_list});
load_bbox_grid_sub = cat(1, grid_info.bbox_grid_sub{layer_list});
num_load_bbox = size(load_bbox_grid_sub, 1);

[wb_link_feature_cell, wb_node_feature_cell,...
    wb_mask_1_proj_cell, wb_mask_2_proj_cell, wb_mask_3_proj_cell, ...
    wb_ai_all_vw, wb_ai_cap_vw, wb_ai_all_lw, wb_ai_cap_lw, ...
    wb_cube_stat] = deal(cell(num_load_bbox, 1));
%%
pool_obj = gcp('nocreate');
delete(pool_obj);
parpool(num_core);
parfor iter_bbox = 1 : num_load_bbox
    tmp_grid_sub = load_bbox_grid_sub(iter_bbox, :);
    % Load data
    try
        tmp_mask_str = DataManager.load_block_mask(dataset_name, stack, recon_mask_name, ...
            tmp_grid_sub(1), tmp_grid_sub(2), tmp_grid_sub(3));
    catch ME
        fprintf('Unable to load reconstructed mask (%d, %d, %d)\n', tmp_grid_sub);
        fprintf('Error message: %s\n', ME.message);
        continue;
    end
    try
        node_features = tmp_mask_str.node.features;
        link_features = tmp_mask_str.link.features;
        tmp_stat_str = tmp_mask_str.stat;
        %% Further processing features here
        node_features = fun_analysis_postprocess_node_features(node_features, true);
        link_features = fun_analysis_postprocess_link_features(link_features);
%         tmp_stat_str.node_density_mm3 = size(node_features, 1) ./ 0.24^3;
        %%
        wb_link_feature_cell{iter_bbox} = link_features;
        wb_node_feature_cell{iter_bbox} = node_features;
        wb_mask_1_proj_cell{iter_bbox} = tmp_mask_str.max_proj_1;
        wb_mask_2_proj_cell{iter_bbox} = tmp_mask_str.max_proj_2;
        wb_mask_3_proj_cell{iter_bbox} = tmp_mask_str.max_proj_3;
        
        if isfield(tmp_stat_str, 'anisotropy_all_vw')
            if isfield(tmp_stat_str.anisotropy_all_vw, ...
                    'fractional_anisotropy')
                wb_ai_all_vw{iter_bbox} = tmp_stat_str.anisotropy_all_vw;
            end
            tmp_stat_str = rmfield(tmp_stat_str, 'anisotropy_all_vw');
        end
        if isfield(tmp_stat_str, 'anisotropy_all_lw')
            if isfield(tmp_stat_str.anisotropy_all_lw, ...
                    'fractional_anisotropy')
                wb_ai_all_lw{iter_bbox} = tmp_stat_str.anisotropy_all_lw;
            end
            tmp_stat_str = rmfield(tmp_stat_str, 'anisotropy_all_lw');
        end
        if isfield(tmp_stat_str, 'anisotropy_capillary_vw')
            if isfield(tmp_stat_str.anisotropy_capillary_vw, ...
                    'fractional_anisotropy')
                wb_ai_cap_vw{iter_bbox} = tmp_stat_str.anisotropy_capillary_vw;
            end
            tmp_stat_str = rmfield(tmp_stat_str, 'anisotropy_capillary_vw');
        end
        if isfield(tmp_stat_str, 'anisotropy_capillary_lw')
            if isfield(tmp_stat_str.anisotropy_capillary_lw, ...
                    'fractional_anisotropy')
                wb_ai_cap_lw{iter_bbox} = tmp_stat_str.anisotropy_capillary_lw;
            end
            tmp_stat_str = rmfield(tmp_stat_str, 'anisotropy_capillary_lw');
        end
        wb_cube_stat{iter_bbox} = tmp_stat_str;
    catch ME
        fprintf('Unable to process cube (%d, %d, %d)\n', tmp_grid_sub);
        rethrow(ME);
    end
end
pool_obj = gcp('nocreate');
delete(pool_obj);
%% Save to structure
wholebrain_stat_str = struct;
wholebrain_stat_str.grid_size = grid_info.grid_size;
wholebrain_stat_str.grid_ind = load_bbox_grid_ind;
wholebrain_stat_str.grid_sub = load_bbox_grid_sub;
wholebrain_stat_str.link_features = wb_link_feature_cell;
wholebrain_stat_str.node_features = wb_node_feature_cell;
wholebrain_stat_str.mask_max_proj_1 = wb_mask_1_proj_cell;
wholebrain_stat_str.mask_max_proj_2 = wb_mask_2_proj_cell;
wholebrain_stat_str.mask_max_proj_3 = wb_mask_3_proj_cell;

wholebrain_stat_str.wb_ai_all_vw = wb_ai_all_vw;
wholebrain_stat_str.wb_ai_all_lw = wb_ai_all_lw;
wholebrain_stat_str.wb_ai_cap_vw = wb_ai_cap_vw;
wholebrain_stat_str.wb_ai_cap_lw = wb_ai_cap_lw;
wholebrain_stat_str.cube_stat = wb_cube_stat;
end