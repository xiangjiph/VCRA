function exit_code = fun_simulation_OT_internal_subgrid(grid_c_info, load_skl_name, grid_c_label, opt)
% fun_analysis_internal_subgrid is a wrap up function for solving
% Poisson equation in the reconstructed vessel mask and analyze the
% dependence of oxygen partial pressure on local vessel network geometry

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
%% Parse input
dataset_name = grid_c_info.dataset_name;
stack = grid_c_info.stack;
grid_info = grid_c_info.grid_ori;

data_folder_name = opt.save_folder_name;
use_gpu_Q = opt.use_gpu_Q;
% Parameters
recon_max_error_rate = opt.recon_max_error_rate;
int_expand_half_length = opt.recon_pad_half_length; % before downsampling
downsample_rate = opt.recon_downsample_rate;
local_max_window_size_ds = opt.local_extrema_window_size ./ downsample_rate;
int_expand_half_length_ds = int_expand_half_length ./ downsample_rate;

if isfield(opt, 'min_recon_vessel_radius_um')
    min_recon_radius = opt.min_recon_vessel_radius_um;
else
    min_recon_radius = 0;
end

inhomogeneous_term = opt.inhomogeneous_term;
krogh_coeff = inhomogeneous_term / 2;
r_cap = opt.krogh_vessel_r_um;
%%
tmp_grid_c_sub = grid_c_info.bbox_grid_sub_list(grid_c_label, :);
tmp_grid_c_ind = sub2ind(grid_c_info.grid_size, tmp_grid_c_sub(1), ...
    tmp_grid_c_sub(2), tmp_grid_c_sub(3));
assert(grid_c_info.bbox_xyz_label_array(tmp_grid_c_ind) == grid_c_label, 'Incorrect grid ind');
if isempty(grid_c_info.internal_subgrid_valid_idx{tmp_grid_c_ind})
    fprintf('No valid internal sub-grid in this combined grid\n');
    exit_code = 1;
    return
end
%% Load all the subgrid data
% All the subgrid inside the combined grid
grid_c_subgrid_sub = grid_c_info.sub_grid_sub{tmp_grid_c_ind};
grid_c_num_subgrid = size(grid_c_subgrid_sub, 1);
% Internal subgrids inside the combined grid
int_subgrid_bbox_mmxx_list = grid_c_info.internal_subgrid_bbox_mmxx{tmp_grid_c_ind};
int_bbox_exp_mm = max([1,1,1], min(int_subgrid_bbox_mmxx_list(:, 1:3), [], 1) - int_expand_half_length);
int_bbox_exp_xx = min(grid_c_info.data_size, max(int_subgrid_bbox_mmxx_list(:, 4:6), [], 1) + int_expand_half_length);
int_bbox_exp_size = int_bbox_exp_xx - int_bbox_exp_mm + 1;
int_bbox_local_bbox = int_subgrid_bbox_mmxx_list - [int_bbox_exp_mm, int_bbox_exp_mm] + 1;

grid_c_subgrid_valid_array = grid_c_info.sub_grid_label_array{tmp_grid_c_ind} > 0;

assert(nnz(grid_c_subgrid_valid_array) == grid_c_num_subgrid, 'Number of valid subgrid does not match');
% Load skeleton and combined them
grid_c_skl_ind = cell(grid_c_num_subgrid, 1);
grid_c_r = cell(grid_c_num_subgrid, 1);
for iter_cube = 1 : grid_c_num_subgrid
    tmp_idx_1 = grid_c_subgrid_sub(iter_cube, 1);
    tmp_idx_2 = grid_c_subgrid_sub(iter_cube, 2);
    tmp_layer = grid_c_subgrid_sub(iter_cube, 3);
    try
        tmp_skel = DataManager.load_block_skl(dataset_name, stack, load_skl_name, tmp_idx_1, tmp_idx_2, tmp_layer);
    catch
        fprintf('Fail to read the skeleton file of cube (%d, %d, %d) . Skip this block\n', tmp_idx_1, tmp_idx_2, tmp_layer);
        continue;
    end
    tmp_skl_sub = fun_ind2sub(tmp_skel.block_size, tmp_skel.ind);
    tmp_skl_sub_g = bsxfun(@plus, tmp_skl_sub, tmp_skel.global_bbox_mmll(1:3) - 1);
    tmp_skl_sub_int_exp = tmp_skl_sub_g - int_bbox_exp_mm + 1;
    
    tmp_in_bbox_Q = all(tmp_skl_sub_int_exp > 0, 2) & ...
        all( bsxfun(@le, tmp_skl_sub_int_exp, int_bbox_exp_size), 2);
    tmp_skl_sub_valid = tmp_skl_sub_int_exp(tmp_in_bbox_Q , : );
    
    grid_c_skl_ind{iter_cube} = sub2ind(int_bbox_exp_size, tmp_skl_sub_valid(:,1), ...
        tmp_skl_sub_valid(:, 2), tmp_skl_sub_valid(:,3));
    grid_c_r{iter_cube} = tmp_skel.r(tmp_in_bbox_Q);
end
grid_c_skl_ind = cat(1, grid_c_skl_ind{:});
grid_c_r = cat(1, grid_c_r{:});
assert(numel(grid_c_skl_ind) == numel(grid_c_r), 'Number of skeleton voxels does not equal to the number of radius');

[grid_c_skl_ind, tmp_unique_idx, ~] = unique(grid_c_skl_ind);
grid_c_r = grid_c_r(tmp_unique_idx);

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
%% Reconstruction
vessel_recon_mask = fun_graph_to_reconstruction_mask(vessel_graph, false, recon_max_error_rate);
system_size_ds = round(size(vessel_recon_mask) / downsample_rate);
%% Solve Poisson equation
pO2_n_result_ds2 = fun_simulation_OT_solve_multigrid_iter(vessel_recon_mask, ...
    downsample_rate, inhomogeneous_term, use_gpu_Q);
%% Local distance transform properties
assert(all(pO2_n_result_ds2.vessel_mask == system_size_ds), 'Mismatched system size');
mask_rz_dt = bwdist(pO2_n_result_ds2.vessel_mask) .* downsample_rate;

dt_max_info = fun_analysis_get_local_extrema_info(mask_rz_dt, local_max_window_size_ds, 'max');

is_valid_dt_lm_Q = fun_voxel_sub_in_bbox_mmxx_Q(dt_max_info.sub, ...
    [repelem(int_expand_half_length_ds, 1, 3) + 1, system_size_ds - int_expand_half_length_ds]) & ...
    ~pO2_n_result_ds2.vessel_mask(dt_max_info.ind);

dt_max_info = fun_structure_field_indexing(dt_max_info, is_valid_dt_lm_Q);
dt_max_info.pO2_v = pO2_n_result_ds2.pO2_array(dt_max_info.ind);
%% Local pO2 properties
pO2_min_info = fun_analysis_get_local_extrema_info(pO2_n_result_ds2.pO2_array, local_max_window_size_ds, 'min');

is_valid_pO2_lm_Q = fun_voxel_sub_in_bbox_mmxx_Q(pO2_min_info.sub, ...
    [repelem(int_expand_half_length_ds, 1, 3) + 1, system_size_ds - int_expand_half_length_ds]) & ...
    ~pO2_n_result_ds2.vessel_mask(pO2_min_info.ind);

pO2_min_info = fun_structure_field_indexing(pO2_min_info, is_valid_pO2_lm_Q);
pO2_min_info.dt_v = mask_rz_dt(pO2_min_info.ind);
%% Analyze each internal grid
int_bbox_local_bbox_ds = ceil(int_bbox_local_bbox ./ downsample_rate);
int_bbox_grid_label = grid_c_info.internal_subgrid_label{tmp_grid_c_ind};
num_int_cube = numel(int_bbox_grid_label);
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
    %% Record numerical solution information
    tmp_cube_str.final_maximum_update = pO2_n_result_ds2.final_maximum_update;
    tmp_cube_str.final_iteration = pO2_n_result_ds2.final_iteration;
    %% Regional statistics - seperate stat
    tmp_cube_str.local_dt_stat = fun_analysis_get_basic_statistics(tmp_int_dt, true);
    tmp_cube_str.local_dt_stat.prtl2val_itp = fun_analysis_get_prtl_to_value_interpolation(tmp_cube_str.local_dt_stat);
    tmp_cube_str.local_dt_stat.val2ptrl_itp = fun_analysis_get_value_to_prtl_interpolation(tmp_cube_str.local_dt_stat);
    
    tmp_cube_str.local_pO2_stat = fun_analysis_get_basic_statistics(tmp_int_pO2);
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
    tmp_int_lm_dt_Q = fun_voxel_sub_in_bbox_mmxx_Q(dt_max_info.sub, tmp_local_bbox_mmxx);
    tmp_cube_str.dt_lm = fun_structure_field_indexing(dt_max_info, tmp_int_lm_dt_Q);
        % Transform the coordinate
    tmp_cube_str.dt_lm.sub = ceil((tmp_cube_str.dt_lm.sub - tmp_local_bbox_mmxx(1:3) + 1) .* downsample_rate); % Scale back
    tmp_cube_str.dt_lm.ind = sub2ind(tmp_cube_str.block_size, tmp_cube_str.dt_lm.sub(:, 1), ...
        tmp_cube_str.dt_lm.sub(:, 2), tmp_cube_str.dt_lm.sub(:, 3));
    %% Local oxygen minimum inside the cube
    tmp_int_lm_pO2_Q = fun_voxel_sub_in_bbox_mmxx_Q(pO2_min_info.sub, tmp_local_bbox_mmxx);
    tmp_cube_str.pO2_lm = fun_structure_field_indexing(pO2_min_info, tmp_int_lm_pO2_Q);
        % Transform the coordinate
    tmp_cube_str.pO2_lm.sub = ceil((tmp_cube_str.pO2_lm.sub - tmp_local_bbox_mmxx(1:3) + 1) .* downsample_rate); 
    tmp_cube_str.pO2_lm.ind = sub2ind(tmp_cube_str.block_size, tmp_cube_str.pO2_lm.sub(:, 1), ...
        tmp_cube_str.pO2_lm.sub(:, 2), tmp_cube_str.pO2_lm.sub(:, 3));
    %% Pair the local maxima of DT with local minimum of pO2
    tmp_pdist_oxy_2_dt = pdist2(tmp_cube_str.pO2_lm.sub, tmp_cube_str.dt_lm.sub) .* downsample_rate;
    [tmp_cube_str.paired_extrema.pO2_list_idx, tmp_cube_str.paired_extrema.dt_list_idx,...
        tmp_cube_str.paired_extrema.dist] = fun_find_col_row_co_minimum(tmp_pdist_oxy_2_dt);
    %% Fit the effective Krogh model against voxel pO2 vs dt
%     tmp_local_dt_lm_med = median(tmp_cube_str.dt_lm.v);
    tmp_local_dt_lm_med = median(tmp_cube_str.pO2_lm.dt_v);
    if isnan(tmp_local_dt_lm_med)
        tmp_local_dt_lm_med = median(tmp_cube_str.dt_lm.v);
        tmp_cube_str.fit_Krogh.d_max_type = 'med_dt_lm_v';
    else
        tmp_cube_str.fit_Krogh.d_max_type = 'med_pO2_dt_v';
    end
    tmp_krogh_fun = @(x) (- krogh_coeff * ((tmp_local_dt_lm_med + r_cap).^ 2 .* log((x + r_cap) ./ r_cap) - ...
            ((x + r_cap).^ 2 - r_cap^2)/2));
    tmp_selected_fit_Q = (tmp_int_dt <= tmp_local_dt_lm_med);
    tmp_cube_str.fit_Krogh.fit_fun_hdl = tmp_krogh_fun;
    if any(tmp_selected_fit_Q)
        lin_fit_hdl = fitlm(tmp_krogh_fun(tmp_int_dt(tmp_selected_fit_Q)), tmp_int_pO2(tmp_selected_fit_Q), 'Intercept', false);
        tmp_cube_str.fit_Krogh.corr_coeff = lin_fit_hdl.Coefficients.Estimate(1);
        tmp_cube_str.fit_Krogh.corr_coeff_SE = lin_fit_hdl.Coefficients.SE(1);
        tmp_cube_str.fit_Krogh.Rsquared = lin_fit_hdl.Rsquared;
    else
        [tmp_cube_str.fit_Krogh.corr_coeff, tmp_cube_str.fit_Krogh.corr_coeff_SE, ...
            tmp_cube_str.fit_Krogh.Rsquared] = deal(nan);
    end
    tmp_cube_str.fit_Krogh.d_max = tmp_local_dt_lm_med;
    tmp_cube_str.fit_Krogh.cap_r = r_cap;
    %% Fit the effective Krogh model against median / mean curve? 
    %% Save result
    DataManager.write_analysis_data_in_grid(tmp_cube_str, data_folder_name, ...
        dataset_name, stack, tmp_cube_str.grid_version, tmp_cube_str.grid_label);
end
exit_code = 0;
end