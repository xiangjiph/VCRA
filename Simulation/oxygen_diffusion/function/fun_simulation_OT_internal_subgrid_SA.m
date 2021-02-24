function [exit_code, varargout] = fun_simulation_OT_internal_subgrid_SA(grid_c_info, load_skl_name, grid_c_label, opt)

%% Parameters
if ~isa(opt, 'struct')
    if isfile(opt)
        opt = load(opt);
    else
        error('The input opt should be either a valid path to the mat file or a MATLAB structure');
    end
end
persistent DataManager
if isempty(DataManager)
    DataManager = FileManager;
end
if nargout == 2
    output_cube_stat_Q = true;
    varargout = [];
else
    output_cube_stat_Q = false;
end
%% Parse input
dataset_name = grid_c_info.dataset_name;
stack = grid_c_info.stack;
grid_info = grid_c_info.grid_ori;

data_folder_name = opt.save_folder_name;
use_gpu_Q = opt.use_gpu_Q;
mask_itp_method = opt.mask_itp_method;

recon_max_error_rate = opt.recon_max_error_rate;
int_expand_half_length = opt.recon_pad_half_length; % before downsampling
lm_wz_list_um = opt.local_extrema_window_size; % before downsampling

if isfield(opt, 'min_recon_vessel_radius_um')
    min_recon_radius = opt.min_recon_vessel_radius_um;
else
    min_recon_radius = 0;
end

inhomogeneous_term = opt.inhomogeneous_term;
downsample_rate = opt.recon_downsample_rate;
lm_wz_list_pxl = lm_wz_list_um ./ downsample_rate;
int_expand_half_length_ds = int_expand_half_length ./ downsample_rate;
%%
grid_c_sub = grid_c_info.bbox_grid_sub_list(grid_c_label, :);
grid_c_ind = sub2ind(grid_c_info.grid_size, grid_c_sub(1), ...
    grid_c_sub(2), grid_c_sub(3));
assert(grid_c_info.bbox_xyz_label_array(grid_c_ind) == grid_c_label, 'Incorrect grid ind');
if isempty(grid_c_info.internal_subgrid_valid_idx{grid_c_ind})
    fprintf('No valid internal sub-grid in this combined grid\n');
    exit_code = 1;
    return
end
%% Load all the subgrid data
tic_load_file = tic;
[grid_c_skl_ind, grid_c_r, int_bbox_local_bbox_mmxx, int_bbox_exp_size] = ...
    fun_simulation_OT_load_skeleton_in_grid(grid_c_info, load_skl_name, grid_c_ind, int_expand_half_length);
fprintf('Finish loading skeleton. Elapsed time is %f seconds.\n', toc(tic_load_file));
if isempty(grid_c_skl_ind)
   fprintf('Empty block! Return\n');
   exit_code = 2;
   return;
end
%% Convert to graph
% Modify vessel radius for reconstruction 
grid_c_r = max(grid_c_r, min_recon_radius);
vessel_graph = fun_skeleton_to_graph(grid_c_skl_ind, int_bbox_exp_size);
vessel_graph.radius = sparse(grid_c_skl_ind, ones(vessel_graph.num.skeleton_voxel, 1), ...
    double(grid_c_r), vessel_graph.num.block_voxel, 1);
%% Reconstruction vessel mask at 1 um resolution
vessel_recon_mask = fun_graph_to_reconstruction_mask(vessel_graph, false, recon_max_error_rate);
%% Solve Poisson equation
pO2_n_result_ds2 = fun_simulation_OT_solve_multigrid_iter(vessel_recon_mask, ...
    downsample_rate, inhomogeneous_term, mask_itp_method, use_gpu_Q);
%% Local distance transform properties
% local maxima found by larger window size are a subset of local maxima
% found by smaller window size
mask_size_ds = size(pO2_n_result_ds2.vessel_mask);
valid_mask_bbox_mmxx = [repelem(int_expand_half_length_ds, 1, 3) + 1, ...
    mask_size_ds - int_expand_half_length_ds];
valid_extrema_mask = false(mask_size_ds);
valid_extrema_mask(valid_mask_bbox_mmxx(1) : valid_mask_bbox_mmxx(4), ...
    valid_mask_bbox_mmxx(2) : valid_mask_bbox_mmxx(5), ...
    valid_mask_bbox_mmxx(3) : valid_mask_bbox_mmxx(6)) = true;
valid_extrema_mask = valid_extrema_mask & ~pO2_n_result_ds2.vessel_mask;

tic_search_dtlm = tic;
mask_rz_dt = bwdist(pO2_n_result_ds2.vessel_mask) .* downsample_rate;

[dt_max_info, pO2_min_info] = fun_simulation_OT_SA_wdw_sz(mask_rz_dt, pO2_n_result_ds2.pO2_array, ...
    lm_wz_list_pxl, valid_extrema_mask);
fprintf('Finish searching for local extrema. Elapsed time is %f secodns.\n', toc(tic_search_dtlm));
%% Analyze each internal grid
int_bbox_local_bbox_ds = ceil(int_bbox_local_bbox_mmxx ./ downsample_rate);
int_bbox_grid_label = grid_c_info.internal_subgrid_label{grid_c_ind};
num_int_cube = numel(int_bbox_grid_label);
if output_cube_stat_Q
    cube_stat_cell = cell(num_int_cube, 1);
end
tic_write = tic;
for iter_int_cube = 1 : num_int_cube
    tmp_cube_label = int_bbox_grid_label(iter_int_cube);
    tmp_cube_str = fun_grid_get_single_cube_info(grid_info, tmp_cube_label);
    tmp_local_bbox_mmxx = int_bbox_local_bbox_ds(iter_int_cube, :);
    tmp_local_bbox_mmll = [tmp_local_bbox_mmxx(1:3), tmp_local_bbox_mmxx(4:6) - tmp_local_bbox_mmxx(1:3) + 1];
    % Crop the local dt and pO2 array
    tmp_int_dt = crop_bbox3(mask_rz_dt, tmp_local_bbox_mmll);
    tmp_int_pO2 = crop_bbox3(pO2_n_result_ds2.pO2_array, tmp_local_bbox_mmll);
    % Select the voxel outside the vessel mask 
    tmp_int_tissue_mask = (tmp_int_dt ~= 0);
    tmp_int_dt = tmp_int_dt(tmp_int_tissue_mask);
    tmp_int_pO2 = tmp_int_pO2(tmp_int_tissue_mask);
    %% Regional statistics - seperate stat
    tmp_cube_str.local_dt_stat = fun_analysis_get_basic_statistics(tmp_int_dt, true);
    tmp_cube_str.local_dt_stat.prtl2val_itp = fun_analysis_get_prtl_to_value_interpolation(tmp_cube_str.local_dt_stat);
    tmp_cube_str.local_dt_stat.val2ptrl_itp = fun_analysis_get_value_to_prtl_interpolation(tmp_cube_str.local_dt_stat);
    
    tmp_cube_str.local_pO2_stat = fun_analysis_get_basic_statistics(tmp_int_pO2, true);    
    tmp_cube_str.local_pO2_stat.prtl2val_itp = fun_analysis_get_prtl_to_value_interpolation(tmp_cube_str.local_pO2_stat);
    tmp_cube_str.local_pO2_stat.val2ptrl_itp = fun_analysis_get_value_to_prtl_interpolation(tmp_cube_str.local_pO2_stat);
    %% Regional statistics - joint stat
    tmp_dt_edge = 0.5 : tmp_cube_str.local_dt_stat.max;
    tmp_cube_str.pO2_stat_in_dt_bin = fun_analysis_get_y_stat_in_x_bin(tmp_int_dt, ...
        tmp_int_pO2, tmp_dt_edge);
    %% 2D histogram - count - normalization can be done later if needed
    [tmp_cube_str.pO2_dt_hist2.count, tmp_cube_str.pO2_dt_hist2.dt_edge, ...
        tmp_cube_str.pO2_dt_hist2.pO2_edge] = histcounts2(tmp_int_dt, tmp_int_pO2);
    %% Local DT maxima inside the cube
    tmp_cube_str.local_extrema_window_size = lm_wz_list_um;
    tmp_cube_str.dt_lm = fun_simulation_OT_SA_get_extrema_in_bbox(dt_max_info, ...
        tmp_local_bbox_mmxx, downsample_rate);
    [tmp_cube_str.dt_lm_stat.dt_mean, tmp_cube_str.dt_lm_stat.dt_median] = ...
        fun_analysis_get_mean_N_median_in_each_cell(tmp_cube_str.dt_lm.dt);
    [tmp_cube_str.dt_lm_stat.pO2_mean, tmp_cube_str.dt_lm_stat.pO2_median] = ...
        fun_analysis_get_mean_N_median_in_each_cell(tmp_cube_str.dt_lm.pO2);
    %% Local oxygen minimum inside the cube
    tmp_cube_str.pO2_lm = fun_simulation_OT_SA_get_extrema_in_bbox(pO2_min_info, ...
        tmp_local_bbox_mmxx, downsample_rate);
    [tmp_cube_str.pO2_lm_stat.dt_mean, tmp_cube_str.pO2_lm_stat.dt_median] = ...
        fun_analysis_get_mean_N_median_in_each_cell(tmp_cube_str.pO2_lm.dt);
    [tmp_cube_str.pO2_lm_stat.pO2_mean, tmp_cube_str.pO2_lm_stat.pO2_median] = ...
        fun_analysis_get_mean_N_median_in_each_cell(tmp_cube_str.pO2_lm.pO2);
    % Add the percentile of the dt_lm and pO2_lm    
    %% Save result
    DataManager.write_analysis_data_in_grid(tmp_cube_str, data_folder_name, ...
        dataset_name, stack, tmp_cube_str.grid_version, tmp_cube_str.grid_label);
    if output_cube_stat_Q
       cube_stat_cell{iter_int_cube} = tmp_cube_str; 
    end
end
fprintf('Finish writing files. Elapsed time is %f seconds.\n', toc(tic_write));
exit_code = 0;
if output_cube_stat_Q
    varargout{1} = cube_stat_cell;
end
end