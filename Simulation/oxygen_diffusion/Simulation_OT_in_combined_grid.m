set_env;
dataset_name = 'WholeBrain';
stack = 'ML_2018_08_15';
load_skl_name = '240_cube_rec';
combined_grid_name = '240_cube_combined_5_o_2';
grid_c_info = DataManager.load_grid(dataset_name, stack, combined_grid_name);
grid_info = grid_c_info.grid_ori;
output_graph_name = 'OT_simulation';
data_folder_name = 'pO2';
overwrite_Q = false;
vis_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), 'Tissue_pO2_simulation');
% Might be worth to apply a minimum capillary radius for: 
% 1. Downsampling accuracy
% 2. Reflect the real capillary radius in vivo
min_cap_r = 2; % um 
%% 
recon_max_error_rate = 0.1;
int_expand_half_length = 80; % before downsampling
local_max_window_size = 40; % before downsampling
downsample_rate = 2;
local_max_window_size_ds = local_max_window_size ./ downsample_rate;
int_expand_half_length_ds = int_expand_half_length ./ downsample_rate;

r_cap = 2;
krogh_coeff = 1/2;
est_corr_coeff = 0.5;
gpuDevice_label = 1;
gpuDevice(gpuDevice_label);
inhomogeneous_term = 1;
%% Test input
grid_c_label = 860;
%%
[grid_c_skl_ind, grid_c_r, int_bbox_local_bbox_mmxx, int_bbox_exp_size] = ...
    fun_simulation_OT_load_skeleton_in_grid(grid_c_info, load_skl_name, grid_c_label, int_expand_half_length);
%% Convert to graph
vessel_graph = fun_skeleton_to_graph(grid_c_skl_ind, int_bbox_exp_size);
vessel_graph.radius = sparse(grid_c_skl_ind, ones(vessel_graph.num.skeleton_voxel, 1), ...
    max(min_cap_r, double(grid_c_r)), vessel_graph.num.block_voxel, 1);
%% Reconstruction
vessel_recon_mask = fun_graph_to_reconstruction_mask(vessel_graph, false, recon_max_error_rate);
%% Downsample the mask
int_bbox_local_bbox_ds = ceil(int_bbox_local_bbox_mmxx ./ downsample_rate);
tmp_target_size = round(int_bbox_exp_size ./ downsample_rate);
tmp_recon_rz = imresize3(uint8(vessel_recon_mask), tmp_target_size) > 0;
mask_rz_dt = bwdist(tmp_recon_rz) .* downsample_rate;
%% Local distance transform properties
% local maxima found by larger window size are a subset of local maxima
% found by smaller window size
dt_max_info = fun_analysis_get_local_extrema_info(mask_rz_dt, local_max_window_size_ds, 'max');
is_valid_dt_lm_Q = fun_voxel_sub_in_bbox_mmxx_Q(dt_max_info.sub, ...
    [repelem(int_expand_half_length_ds, 1, 3) + 1, tmp_target_size - int_expand_half_length_ds]) & ...
    ~tmp_recon_rz(dt_max_info.ind);

dt_max_info = fun_structure_field_indexing(dt_max_info, is_valid_dt_lm_Q & is_not_zero_Q);
%% Estimate pO2 initial condition
dt_lm_v_med = median(dt_max_info.v);
krogh_fun = @(x) (- krogh_coeff * ((dt_lm_v_med + r_cap).^ 2 .* log((x + r_cap) ./ r_cap) - ...
            ((x + r_cap).^ 2 - r_cap^2)/2));
tmp_ini_pO2_n = est_corr_coeff * krogh_fun(min(mask_rz_dt, dt_lm_v_med));
%% Solve Poisson equation
tmp_inhomgeneous_term = inhomogeneous_term * downsample_rate.^ 2;
pO2_n_result_ds2 = fun_simulation_OT_solve_ct_diff_itr_n(tmp_recon_rz, tmp_inhomgeneous_term, tmp_ini_pO2_n, true);    
%% Local pO2 properties
pO2_min_info = fun_analysis_get_local_extrema_info(pO2_n_result_ds2.pO2_array, local_max_window_size_ds, 'min');

is_valid_pO2_lm_Q = fun_voxel_sub_in_bbox_mmxx_Q(pO2_min_info.sub, ...
    [repelem(int_expand_half_length_ds, 1, 3) + 1, tmp_target_size - int_expand_half_length_ds]) & ...
    ~tmp_recon_rz(pO2_min_info.ind);

pO2_min_info = fun_structure_field_indexing(pO2_min_info, is_valid_pO2_lm_Q);
pO2_min_info.dt_v = mask_rz_dt(pO2_min_info.ind);
dt_max_info.pO2_v = pO2_n_result_ds2.pO2_array(dt_max_info.ind);
%% Analyze each internal grid
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
    tmp_pdist_oxy_2_dt = pdist2(tmp_cube_str.pO2_lm.sub, tmp_cube_str.dt_lm.sub);
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
    %% Save result
    exit_code = DataManager.write_analysis_data_in_grid(tmp_cube_str, ...
        data_folder_name, dataset_name, stack, tmp_cube_str.grid_version, ...
        tmp_cube_str.grid_label);
end
%% Visualization
tmp_int_dt = crop_bbox3(mask_rz_dt, tmp_local_bbox_mmll);
tmp_int_pO2 = crop_bbox3(pO2_n_result_ds2.pO2_array, tmp_local_bbox_mmll);
% tmp_int_dt = mask_rz_dt;
% tmp_int_pO2 = pO2_n_result_ds2.pO2_array;
%% Video: pO2 and distance transform
video_fp = fullfile(DataManager.fp_visualization_folder(dataset_name, ...
    stack), 'paper', sprintf('%s_%s_pO2_dt_vis_%d.mp4', dataset_name, stack, ...
    tmp_cube_str.grid_ind));
video_folder = fileparts(video_fp);
if ~isfolder(video_folder)
    mkdir(video_folder);
end

avi_str = VideoWriter(video_fp);
avi_str.FrameRate = 5;
avi_str.Quality = 90;
open(avi_str);
fig_hdl = figure('Visible', 'on');
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [2, 1];
ax_hdl_2 = subplot(1,2,1);
ax_hdl_2 = subplot(1,2,2);
for vis_sec = 1 : size(tmp_int_pO2, 3)
    imagesc(ax_hdl_2, tmp_int_pO2(:, :, vis_sec));
    ax_hdl_2.DataAspectRatio = [1,1,1];
    ax_hdl_2.FontSize = 14;
    ax_hdl_2.XTickLabel = arrayfun(@(x) num2str(x * 2, '%d'), ax_hdl_2.XTick, 'UniformOutput', false);
    ax_hdl_2.YTickLabel = arrayfun(@(x) num2str(x * 2, '%d'), ax_hdl_2.YTick, 'UniformOutput', false);
    ax_hdl_2.XLabel.String = 'X (\mum)';
    ax_hdl_2.YLabel.String = 'Y (\mum)';
    ax_hdl_2.CLim = [-800, 0];
    cbar_hdl_1 = colorbar(ax_hdl_2);
    cbar_hdl_1.Label.String = 'Normalized pO_2';

    ax_hdl_2.Title.String = sprintf('Normalized pO_2 (Z = %d \\mum)', vis_sec * downsample_rate);
    
    tmp_on_plane_pO2_lm_Q = tmp_cube_str.pO2_lm.sub(:, 3) == (vis_sec * downsample_rate);
%     tmp_on_plane_pO2_lm_Q = pO2_min_info.sub(:, 3) == vis_sec;
    if any(tmp_on_plane_pO2_lm_Q)
        hold(ax_hdl_2, 'on')
        tmp_on_plane_sub = tmp_cube_str.pO2_lm.sub(tmp_on_plane_pO2_lm_Q, 1:2) ./ downsample_rate;
%         tmp_on_plane_sub = pO2_min_info.sub(tmp_on_plane_pO2_lm_Q, 1:2);
        scatter(ax_hdl_2, tmp_on_plane_sub(:, 2), tmp_on_plane_sub(:, 1), 'r*');
        hold(ax_hdl_2, 'off');
    end
    
    imagesc(ax_hdl_2, tmp_int_dt(:, :, vis_sec));
    ax_hdl_2.DataAspectRatio = [1,1,1];
    ax_hdl_2.FontSize = 14;
    ax_hdl_2.XTickLabel = arrayfun(@(x) num2str(x * 2, '%d'), ax_hdl_2.XTick, 'UniformOutput', false);
    ax_hdl_2.YTickLabel = arrayfun(@(x) num2str(x * 2, '%d'), ax_hdl_2.YTick, 'UniformOutput', false);
    ax_hdl_2.XLabel.String = 'X (\mum)';
    ax_hdl_2.YLabel.String = 'Y (\mum)';
    cbar_hdl_2 = colorbar(ax_hdl_2);
    ax_hdl_2.CLim = [0, 35];
    cbar_hdl_2.Label.String = 'Distance to the nearest vessel (\mum)';

    ax_hdl_2.Title.String = sprintf('Distance transform (Z = %d \\mum)', vis_sec * downsample_rate);
    tmp_on_plane_dt_lm_Q = tmp_cube_str.dt_lm.sub(:, 3) == (vis_sec * downsample_rate);
%     tmp_on_plane_dt_lm_Q = dt_max_info.sub(:, 3) == vis_sec;
    if any(tmp_on_plane_dt_lm_Q)
        hold(ax_hdl_2, 'on')
        tmp_on_plane_sub = tmp_cube_str.dt_lm.sub(tmp_on_plane_dt_lm_Q, 1:2)  ./ downsample_rate;
%         tmp_on_plane_sub = dt_max_info.sub(tmp_on_plane_dt_lm_Q, 1:2);
        scatter(ax_hdl_2, tmp_on_plane_sub(:, 2), tmp_on_plane_sub(:, 1), 'r+');
        hold(ax_hdl_2, 'off');
    end
    writeVideo(avi_str, getframe(fig_hdl));
end
close(avi_str)
%% Video: pO2, local minimum of pO2, and local maximum of DT
video_fp = fullfile(DataManager.fp_visualization_folder(dataset_name, ...
    stack), 'paper', sprintf('%s_%s_pO2_vis_%d.mp4', dataset_name, stack, ...
    tmp_cube_str.grid_ind));
video_folder = fileparts(video_fp);
if ~isfolder(video_folder)
    mkdir(video_folder);
end

avi_str = VideoWriter(video_fp);
avi_str.FrameRate = 5;
avi_str.Quality = 90;
open(avi_str);
fig_hdl = figure('Visible', 'on');
ax_hdl_2 = axes(fig_hdl);
for vis_sec = 1 : size(tmp_int_pO2, 3)
    imagesc(ax_hdl_2, tmp_int_pO2(:, :, vis_sec));
    ax_hdl_2.DataAspectRatio = [1,1,1];
    ax_hdl_2.FontSize = 14;
    ax_hdl_2.XTickLabel = arrayfun(@(x) num2str(x * 2, '%d'), ax_hdl_2.XTick, 'UniformOutput', false);
    ax_hdl_2.YTickLabel = arrayfun(@(x) num2str(x * 2, '%d'), ax_hdl_2.YTick, 'UniformOutput', false);
    ax_hdl_2.XLabel.String = 'X (\mum)';
    ax_hdl_2.YLabel.String = 'Y (\mum)';
    ax_hdl_2.CLim = [-800, 0];
    cbar_hdl_1 = colorbar(ax_hdl_2);
    cbar_hdl_1.Label.String = 'Normalized pO_2';

    ax_hdl_2.Title.String = sprintf('Normalized pO_2 (Z = %d \\mum)', vis_sec * downsample_rate);
    
    tmp_on_plane_pO2_lm_Q = tmp_cube_str.pO2_lm.sub(:, 3) == (vis_sec * downsample_rate);
%     tmp_on_plane_pO2_lm_Q = pO2_min_info.sub(:, 3) == vis_sec;
    if any(tmp_on_plane_pO2_lm_Q)
        hold(ax_hdl_2, 'on')
        tmp_on_plane_sub = tmp_cube_str.pO2_lm.sub(tmp_on_plane_pO2_lm_Q, 1:2) ./ downsample_rate;
%         tmp_on_plane_sub = pO2_min_info.sub(tmp_on_plane_pO2_lm_Q, 1:2);
        scatter(ax_hdl_2, tmp_on_plane_sub(:, 2), tmp_on_plane_sub(:, 1), 'r*');
        hold(ax_hdl_2, 'off');
    end
    
     tmp_on_plane_dt_lm_Q = tmp_cube_str.dt_lm.sub(:, 3) == (vis_sec * downsample_rate);
%     tmp_on_plane_dt_lm_Q = dt_max_info.sub(:, 3) == vis_sec;
    if any(tmp_on_plane_dt_lm_Q)
        hold(ax_hdl_2, 'on')
        tmp_on_plane_sub = tmp_cube_str.dt_lm.sub(tmp_on_plane_dt_lm_Q, 1:2)  ./ downsample_rate;
%         tmp_on_plane_sub = dt_max_info.sub(tmp_on_plane_dt_lm_Q, 1:2);
        scatter(ax_hdl_2, tmp_on_plane_sub(:, 2), tmp_on_plane_sub(:, 1), 'r+');
        hold(ax_hdl_2, 'off');
    end
    writeVideo(avi_str, getframe(fig_hdl));
    
    if any(tmp_on_plane_pO2_lm_Q) && any(tmp_on_plane_dt_lm_Q)
        tmp_fp = fullfile(DataManager.fp_visualization_folder(dataset_name, ...
            stack), 'paper', sprintf('%s_%s_pO2_vis_%d_sec_%d.png', dataset_name, stack, ...
            tmp_cube_str.grid_ind, vis_sec * 2));
        fun_print_image_in_several_formats(fig_hdl, tmp_fp);
    end
end
close(avi_str)
%% Compare with the direct simulation at the center
% iter_cube = (125 + 1)/2;
% tmp_idx_1 = grid_c_subgrid_sub(iter_cube, 1);
% tmp_idx_2 = grid_c_subgrid_sub(iter_cube, 2);
% tmp_layer = grid_c_subgrid_sub(iter_cube, 3);
% 
% mask_version = '240_cube_recon';
% tmp_recon_str = DataManager.load_block_mask(dataset_name, stack, ...
%     mask_version, tmp_idx_1, tmp_idx_2, tmp_layer);
% subgrid_in_bbox_mmxx = tmp_recon_str.global_bbox_mmxx - [int_bbox_exp_mm, ...
%     int_bbox_exp_mm] + 1;
% subgrid_in_bbox_mm = (subgrid_in_bbox_mmxx(1:3) + 1) ./ 2;
% pO2_grid_c_crop = crop_bbox3(pO2_n_result_ds2.pO2_array, [subgrid_in_bbox_mm, 120, 120, 120]);
% 
% tmp_recon = fun_reconstruct_block_mask(tmp_recon_str);
% tmp_recon_rz_sub = imresize3(uint8(tmp_recon), round(size(tmp_recon) / downsample_rate)) > 0;
% 
% pO2_n_result_ds2_sub = fun_simulation_OT_solve_ct_diff_itr_n(tmp_recon_rz_sub, tmp_inhomgeneous_term, -500, true);  
% pO2_subgrid = pO2_n_result_ds2_sub.pO2_array;
% Visualization
% tmp_sec = 50;
% fig_hdl = figure;
% fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3, 1];
% ax_hdl_1 = subplot(1,3,1);
% imagesc(pO2_grid_c_crop(:, :, tmp_sec));
% ax_hdl_1.Title.String = 'With padding';
% ax_hdl_1.DataAspectRatio = [1,1,1];
% cbar_1 = colorbar(ax_hdl_1);
% ax_hdl_2 = subplot(1,3,2);
% imagesc(pO2_subgrid(:, :, tmp_sec));
% ax_hdl_2.Title.String = 'Without padding';
% ax_hdl_2.DataAspectRatio = [1,1,1];
% cbar_2 = colorbar(ax_hdl_2);
% ax_hdl_2.CLim = ax_hdl_1.CLim;
% ax_hdl_3 = subplot(1,3,3);
% imagesc(pO2_subgrid(:, :, tmp_sec) - pO2_grid_c_crop(:, :, tmp_sec));
% ax_hdl_3.DataAspectRatio = [1,1,1];
% ax_hdl_3.Title.String = 'Without padding - with padding';
% colorbar(ax_hdl_3);
% fig_fp = fullfile(vis_folder, sprintf('%s_%s_boundary_effect_on_pO2_5.png', ...
%     dataset_name, stack));
% fun_print_image_in_several_formats(fig_hdl, fig_fp);
%%
fun_vis_OT_numerical_vs_krogh(tmp_cube_str, true);
    











