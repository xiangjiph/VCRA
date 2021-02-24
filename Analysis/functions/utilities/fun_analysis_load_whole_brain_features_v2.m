function wholebrain_stat_str = fun_analysis_load_whole_brain_features_v2(grid_info, recon_mask_name, layer_list)


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

computed_stat_name = {'mask_volume_ratio', 'mask_surface_area', 'total_length_um', ...
        'total_surface_area_um2', 'total_volume_um3', 'total_capillary_length_um', 'total_capillary_surface_area_um2', ...
        'total_capillary_volume_um3'};
num_computed_stat = numel(computed_stat_name);
anisotropy_field_name = {'anisotropy_all_vw', 'anisotropy_capillary_vw', ...
    'anisotropy_capillary', 'anisotropy_all'};
num_anisotropy_field = numel(anisotropy_field_name);
anisotropy_subfield_name = {'fractional_anisotropy', 'min2max_z', 'svd_1_z'};
num_anisotropy_subfield = numel(anisotropy_subfield_name);

scratch_folder = fullfile(DataManager.Scratch_Folder_Path, 'layer_stat');
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
   
    layer_data = fun_initialized_structure_array_with_fieldname_list(computed_stat_name, sprintf('nan(%d, %d)', grid_info.grid_size(1:2)));
    [layer_data.mask_max_proj_1, layer_data.mask_max_proj_2, layer_data.mask_max_proj_3, ...
        layer_data.link_features, layer_data.node_features] = deal(cell(grid_info.grid_size(1:2)));
    for iter_1 = 1 : num_anisotropy_field
        for iter_2 = 1 : num_anisotropy_subfield
            layer_data.(sprintf('%s_%s', anisotropy_field_name{iter_1}, anisotropy_subfield_name{iter_2})) = nan(grid_info.grid_size(1:2));
        end
    end
    % Load graph in single plane into cell arrays
    for iter_block = 1 : num_valid_graph
        tmp_idx_1 = block_grid_sub(iter_block, 1);
        tmp_idx_2 = block_grid_sub(iter_block, 2);
        tmp_layer = block_grid_sub(iter_block, 3);
        try
            tmp_mask_str = DataManager.load_block_mask(dataset_name, stack, recon_mask_name, ...
                tmp_idx_1, tmp_idx_2, tmp_layer);
            
            tmp_mask_str.block_volume_mm3 = prod(tmp_mask_str.block_size) / 1e9;
            for iter_field_name = 1 : num_computed_stat
                tmp_field_name = computed_stat_name{iter_field_name};
                layer_data.(tmp_field_name)(tmp_idx_1, tmp_idx_2) = tmp_mask_str.(tmp_field_name);
            end
            % Anisotropy
            for iter_1 = 1 : num_anisotropy_field
                if ~isfield(tmp_mask_str.(anisotropy_field_name{iter_1}), anisotropy_subfield_name{1})
                    fprintf('This block (%d, %d, %d) does not have anisotropy data\n', tmp_idx_1, tmp_idx_2, tmp_layer);
                    break;
                else
                    for iter_2 = 1 : num_anisotropy_subfield
                        layer_data.(sprintf('%s_%s', anisotropy_field_name{iter_1}, anisotropy_subfield_name{iter_2}))(tmp_idx_1, tmp_idx_2)...
                            = tmp_mask_str.(anisotropy_field_name{iter_1}).(anisotropy_subfield_name{iter_2});
                    end
                end
            end            
            layer_data.mask_max_proj_1{tmp_idx_1, tmp_idx_2} = tmp_mask_str.max_proj_1;
            layer_data.mask_max_proj_2{tmp_idx_1, tmp_idx_2} = tmp_mask_str.max_proj_2;
            layer_data.mask_max_proj_3{tmp_idx_1, tmp_idx_2} = tmp_mask_str.max_proj_3;
            layer_data.link_features{tmp_idx_1, tmp_idx_2} = tmp_mask_str.link.features;
            layer_data.node_features{tmp_idx_1, tmp_idx_2} = tmp_mask_str.node.features;
        catch ME
            fprintf('Unable to load reconstructed mask (%d, %d, %d)\n', tmp_idx_1, tmp_idx_2, tmp_layer);
            fprintf('Error message: %s\n', ME.message);
            continue;
        end
    end
    % Save to folder
    tmp_save_name = fullfile(scratch_folder, sprintf('%s_%s_layer_%d.mat', dataset_name, stack, grid_layer));
    DataManager.write_data(tmp_save_name, layer_data);
    fprintf('Finish processing layer %d. Elapsed time is %f seconds\n', grid_layer, toc(tmp_tic));
end
%% Load from the scratch folder 
wholebrain_stat_str = fun_initialized_structure_array_with_fieldname_list(computed_stat_name, sprintf('nan(%d, %d, %d)', grid_info.grid_size));
for iter_1 = 1 : num_anisotropy_field
    for iter_2 = 1 : num_anisotropy_subfield
        wholebrain_stat_str.(sprintf('%s_%s', anisotropy_field_name{iter_1}, anisotropy_subfield_name{iter_2})) = nan(grid_info.grid_size);
    end
end
wholebrain_stat_str.grid_size = wb_grid_size;
[wholebrain_stat_str.link_features, wholebrain_stat_str.node_features, ...
    wholebrain_stat_str.mask_max_proj_1, wholebrain_stat_str.mask_max_proj_2, ...
    wholebrain_stat_str.mask_max_proj_3] = deal(cell(wb_grid_size));

for iter_layer = 1 : num_layer_to_process
    tmp_tic = tic;
    grid_layer = layer_list(iter_layer);
    tmp_read_name = fullfile(scratch_folder, sprintf('%s_%s_layer_%d.mat', dataset_name, stack, grid_layer));
    tmp_layer_data = DataManager.load_data(tmp_read_name);
    tmp_field_name = fieldnames(tmp_layer_data);
    for iter_field = 1 : numel(tmp_field_name)
        wholebrain_stat_str.(tmp_field_name{iter_field})(:, :, iter_layer) = tmp_layer_data.(tmp_field_name{iter_field});
    end
    fprintf('Finish loading layer %d. Elapsed time is %f seconds\n', grid_layer, toc(tmp_tic));
end
end