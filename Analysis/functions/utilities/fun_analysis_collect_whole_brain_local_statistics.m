function exit_code = fun_analysis_collect_whole_brain_local_statistics(dataset_name, stack, grid_version, ...
    reconstruction_version, save_folder_name)

persistent DataManager
if isempty(DataManager)
    DataManager = FileManager;
end
grid_info = DataManager.load_grid(dataset_name, stack, grid_version);
%% Options:
save_max_projection_mask_Q = false;
%% Compute 240 cube volume fraction
wb_mask_ds_ratio = 16;
wb_mask = DataManager.load_data(sprintf('%s_mask.nii.gz', fullfile(DataManager.fp_mask_folder(dataset_name, stack, 'whole_brain_d16x_registration'),...
    sprintf('%s_%s_d16x_registration', dataset_name, stack)))) > 0;
is_internal_240_cube_ratio = fun_analysis_get_bbox_in_mask_vol_ratio(wb_mask, grid_info.bbox_xyz_mmxx_list ./ wb_mask_ds_ratio);
grid_info.bbox_volume_ratio_array(grid_info.bbox_grid_ind_list) = is_internal_240_cube_ratio;
%% Load whole brain reconstruction analysis
% Load all the computed features and max projection of the reconstructed images to 3 directions. Takes a few minutes
tic
wholebrain_stat_str = fun_analysis_load_whole_brain_240_cube_features_N_proj(grid_info, reconstruction_version);
fprintf('Finish loading whole brain statistics data. Elapsed time is %f seconds.\n', toc);
%% Compute whole brain statistics
% Takes about 1 minute for 1 layer. Running 18 processes in parallel needs
% about 5 mintues.
tmp_tic = tic;
local_recon_stat = fun_analysis_compute_local_reconstruction_statistics(wholebrain_stat_str);
fprintf('Finish computing local network statistics. Elapsed time is %f seconds.\n', toc(tmp_tic));
% Extract statistics from cell array of structures
scalar_stat_field_name = {'sum', 'mean', 'median', 'num_data', 'std'};
cell_stat_field_name = {'hist_bin_val', 'hist_pdf', 'hist_edge', 'hist_probability'};

[wb_all_vessel, wb_capillary, wb_node] = deal(struct);
for iter_field = 1 : numel(scalar_stat_field_name)
    tmp_field_name = scalar_stat_field_name{iter_field};
    wb_all_vessel.(tmp_field_name) = fun_analysis_get_fields_in_cell_array_str(local_recon_stat.all_link_stat, tmp_field_name);
    wb_capillary.(tmp_field_name) = fun_analysis_get_fields_in_cell_array_str(local_recon_stat.capillary_stat, tmp_field_name);
    wb_node.(tmp_field_name) = fun_analysis_get_fields_in_cell_array_str(local_recon_stat.all_node_stat, tmp_field_name);
end
for iter_field = 1 : numel(cell_stat_field_name)
    tmp_field_name = cell_stat_field_name{iter_field};
    wb_all_vessel.(tmp_field_name) = fun_analysis_get_fields_in_cell_array_str(local_recon_stat.all_link_stat, tmp_field_name, false);
    wb_capillary.(tmp_field_name) = fun_analysis_get_fields_in_cell_array_str(local_recon_stat.capillary_stat, tmp_field_name, false);
    wb_node.(tmp_field_name) = fun_analysis_get_fields_in_cell_array_str(local_recon_stat.all_node_stat, tmp_field_name, false);
end
%% Save max projection in three directions
if save_max_projection_mask_Q
    wb_max_proj_str = struct;
    wb_max_proj_str.dataset_name = dataset_name;
    wb_max_proj_str.stack = stack;
    wb_max_proj_str.grid_version = grid_info.version;
    wb_max_proj_str.grid_info = grid_info;
    wb_max_proj_str.mask_version = reconstruction_version;
    % wb_max_proj_str.filename = sprintf('%s_%s_whole_brian_mask_proj_data.mat', dataset_name, stack);
    % wb_max_proj_str.folder_name = DataManager.fp_analysis_data_folder(dataset_name, stack);
    tmp_cp_field_name_list = {'grid_size', 'grid_ind', 'grid_sub', 'mask_max_proj_1', ...
        'mask_max_proj_2', 'mask_max_proj_3'};
    for iter_field = 1 : numel(tmp_cp_field_name_list)
        tmp_field_name = tmp_cp_field_name_list{iter_field};
        wb_max_proj_str.(tmp_field_name) = wholebrain_stat_str.(tmp_field_name);
    end
    tmp_tic = tic;
    DataManager.write_analysis_data(wb_max_proj_str, dataset_name, stack, sprintf('%s_%s_%s_max_proj_data.mat', ...
        dataset_name, stack, reconstruction_version), save_folder_name);
    fprintf('Finish writing max projection data. Elapsed time is %f seconds.\n', toc(tmp_tic));
end
%% Save 240 cube statistics
% For regional quantification of local link / node feature statistics. 
% Do not need to include the perturbation result as those properties
% depends on segment length. 
wb_240_cube_stat_str = struct;
wb_240_cube_stat_str.dataset_name = dataset_name;
wb_240_cube_stat_str.stack = stack;
wb_240_cube_stat_str.grid_version = grid_info.version;
wb_240_cube_stat_str.mask_version = reconstruction_version;
wb_240_cube_stat_str.cube_in_brain_mask_ratio = is_internal_240_cube_ratio;
wb_240_cube_stat_str.grid_info = grid_info;
tmp_test_str = wholebrain_stat_str.cube_stat{1};
tmp_cp_field_name_list = fieldnames(tmp_test_str);
anisotropy_field_name = {'anisotropy_all_vw', 'anisotropy_all_lw', ...
    'anisotropy_capillary_vw', 'anisotropy_capillary_lw'};
tmp_cp_field_name_list = setdiff(tmp_cp_field_name_list, anisotropy_field_name);
tmp_is_str_Q = ~cellfun(@isempty, wholebrain_stat_str.cube_stat);
for iter_field = 1 : numel(tmp_cp_field_name_list)
    tmp_field_name = tmp_cp_field_name_list{iter_field};
    tmp_data = nan(size(tmp_is_str_Q));
    tmp_data(tmp_is_str_Q) = cellfun(@(x)x.(tmp_field_name), wholebrain_stat_str.cube_stat(tmp_is_str_Q));
    wb_240_cube_stat_str.(tmp_field_name) = tmp_data;
end
% Anisotropy
tmp_ai_extract_field = {'svd_max_vec', 'fractional_anisotropy', 'corr_zr', 'corr_zr_p', ...
    'svd_value_ratio', 'svd_1_z', 'fa_z', 'svd_1_p', 'fa_p'};
wb_ai_field_name = {'wb_ai_all_vw', 'wb_ai_all_lw', 'wb_ai_cap_vw', 'wb_ai_cap_lw'};
for iter_field = 1 : numel(wb_ai_field_name)
    tmp_field_name = wb_ai_field_name{iter_field};
    for iter_ai_field = 1 : numel(tmp_ai_extract_field)
        tmp_ai_field_name = tmp_ai_extract_field{iter_ai_field};
        tmp_data = cell(size(tmp_is_str_Q));
        tmp_elem_data_size = [];
        tmp_is_empty = true(size(tmp_data));
        for iter_bbox = 1 : numel(tmp_data)
            tmp_cube_stat = wholebrain_stat_str.(tmp_field_name){iter_bbox};
            if isfield(tmp_cube_stat, tmp_ai_field_name)
                tmp_elem_data = tmp_cube_stat.(tmp_ai_field_name);
                if iscolumn(tmp_elem_data)
                    tmp_elem_data = tmp_elem_data';
                end
                if isempty(tmp_elem_data_size)
                    tmp_elem_data_size = numel(tmp_elem_data);
                end
                tmp_data{iter_bbox} = tmp_elem_data;
                tmp_is_empty(iter_bbox) = false;
            end
        end
        tmp_array_data = nan(numel(tmp_data), tmp_elem_data_size);
        tmp_array_data(~tmp_is_empty, :) = cat(1, tmp_data{:});
        wb_240_cube_stat_str.(tmp_field_name).(tmp_ai_field_name) = tmp_array_data;
    end
end

% Link, node statistics
cp_link_stat_list = {'length', 'ep2ep_dist', 'dt_median', 'dt_mean', 'dt_cv', 'straightness', 'tortuosity'...
    'surface_area', 'volume', 'is_large_vessel_Q', 'capillary_branching_order', ...
    'shortest_loop_length', 'shortest_loop_geodesic_length', ...
    'nearest_tissue_volume', 'nearest_tissue_dt_mean', 'nearest_tissue_dt_median', ...
    'nearest_tissue_dt_max', 'nearest_tissue_radius', ...
    'noncapillary_nearest_tissue_dt_mean', 'noncapillary_nearest_tissue_dt_max', ...
    'noncapillary_nearest_tissue_dt_median', 'noncapillary_nearest_tissue_dt_volume', ...
    'noncapillary_nearest_capillary_num_vxl', ...
    'dist_to_nearest_noncapillary_mean', 'noncapillary_nearest_tissue_radius', ...
    'dist_to_brain_surface_mean', 'has_no_ep_Q', 'nearest_tissue_2_skl_dist_max'};

cp_node_stat_list = {'nearest_node_dist', 'degree', 'link_length_max', ...
    'link_length_min', 'link_length_median', ...
    'link_ori_ep_sin_elevation', 'link_ep2ep_sin_elevation', 'link_ori_ep_cos_elevation', 'link_ep2ep_cos_elevation', ...
    'path_to_nearest_neighbor_length', 'path_to_nearest_neighbor_geodesic_length', 'dist_to_brain_surface_mean', ...
    'conn_node_drk_min', 'conn_node_drk_med', 'conn_node_drk_max', ...
    'link_ori_agl_min', 'link_ori_agl_med', 'link_ori_agl_max', ...
    'link_ep2ep_agl_min', 'link_ep2ep_agl_med', 'link_ep2ep_agl_max'};
cp_stat_field_list = {'mean', 'median', 'std', 'num_data'};
[wb_240_cube_stat_str.link_all_stat, wb_240_cube_stat_str.link_cap_stat,...
    wb_240_cube_stat_str.node_stat] = deal(struct);
for iter_stat_field = 1 : numel(cp_stat_field_list)
    tmp_stat_field_name = cp_stat_field_list{iter_stat_field};
    for iter_link_feature = 1 : numel(cp_link_stat_list)
        tmp_link_feature_name = cp_link_stat_list{iter_link_feature};
        if isfield(wb_all_vessel.(tmp_stat_field_name), tmp_link_feature_name)
            wb_240_cube_stat_str.link_all_stat.(tmp_stat_field_name).(tmp_link_feature_name) = wb_all_vessel.(tmp_stat_field_name).(tmp_link_feature_name);
        else
            warning('Field does not exist: %s\n', tmp_link_feature_name);
        end
        if isfield(wb_capillary.(tmp_stat_field_name), tmp_link_feature_name)
            wb_240_cube_stat_str.link_cap_stat.(tmp_stat_field_name).(tmp_link_feature_name) = wb_capillary.(tmp_stat_field_name).(tmp_link_feature_name);
        else
            warning('Field does not exist: %s\n', tmp_link_feature_name);
        end
    end
    for iter_node_feautre = 1 : numel(cp_node_stat_list)
        if isfield(wb_node.(tmp_stat_field_name), cp_node_stat_list{iter_node_feautre})
            wb_240_cube_stat_str.node_stat.(tmp_stat_field_name).(cp_node_stat_list{iter_node_feautre}) ...
                = wb_node.(tmp_stat_field_name).(cp_node_stat_list{iter_node_feautre});
        else
            warning('Field does not exist: %s\n', cp_node_stat_list{iter_node_feautre});
        end
    end
end
tmp_tic = tic;
DataManager.write_analysis_data(wb_240_cube_stat_str, dataset_name, stack, sprintf('%s_%s_%s_240_cube_stat_data.mat', ...
    dataset_name, stack, reconstruction_version), save_folder_name);
fprintf('Finish writing 240 cube statistical data. Elpased time is %f seconds.\n', ...
    toc(tmp_tic));
%% Save 240 cube pdf
wb_240_cube_stat_pdf_str = struct;
wb_240_cube_stat_pdf_str.dataset_name = dataset_name;
wb_240_cube_stat_pdf_str.stack = stack;
wb_240_cube_stat_pdf_str.grid_version = grid_info.version;
wb_240_cube_stat_pdf_str.grid_info = grid_info;
wb_240_cube_stat_pdf_str.mask_version = reconstruction_version;
wb_240_cube_stat_pdf_str.cube_in_brain_mask_ratio = is_internal_240_cube_ratio;

cp_stat_field_list = cell_stat_field_name;
for iter_stat_field = 1 : numel(cp_stat_field_list)
    tmp_stat_field_name = cp_stat_field_list{iter_stat_field};
    for iter_link_feature = 1 : numel(cp_link_stat_list)
        tmp_link_feature_name = cp_link_stat_list{iter_link_feature};
        if isfield(wb_all_vessel.(tmp_stat_field_name), tmp_link_feature_name)
            wb_240_cube_stat_pdf_str.link_all_stat.(tmp_stat_field_name).(tmp_link_feature_name) = wb_all_vessel.(tmp_stat_field_name).(tmp_link_feature_name);
        end
        if isfield(wb_capillary.(tmp_stat_field_name), tmp_link_feature_name)
            wb_240_cube_stat_pdf_str.link_cap_stat.(tmp_stat_field_name).(tmp_link_feature_name) = wb_capillary.(tmp_stat_field_name).(tmp_link_feature_name);
        end
    end
    for iter_node_feautre = 1 : numel(cp_node_stat_list)
        if isfield(wb_node.(tmp_stat_field_name), cp_node_stat_list{iter_node_feautre})
            wb_240_cube_stat_pdf_str.node_stat.(tmp_stat_field_name).(cp_node_stat_list{iter_node_feautre}) ...
                = wb_node.(tmp_stat_field_name).(cp_node_stat_list{iter_node_feautre});
        end
    end
end
tmp_tic = tic;
DataManager.write_analysis_data(wb_240_cube_stat_pdf_str, dataset_name, stack, sprintf('%s_%s_%s_240_cube_stat_pdf.mat', ...
    dataset_name, stack, reconstruction_version), save_folder_name);
fprintf('Finish writing 240 cube PDF data. Elpased time is %f seconds.\n', ...
    toc(tmp_tic));

exit_code = 0;
end
